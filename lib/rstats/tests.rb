# A bunch of statistics stuff in ruby.
#
# Copyright (c) 2006 Josh Myer <josh@joshisanerd.com>

# Mostly just a programmatic form of notes from Russell Langley's
# "Practical Statistics Simply Explained", Dover.  ISBN 0486227294

#  The book is highly recommended as a primer, though it is a little
#  light on the theoretical underpinnings. 


require 'rstats'

module RStats
  module Tests
    # The zM test. Langley pp152-9
    #
    # Purpose: The z Test for Measurements (zM) compares a random
    # sample of 1 or more measurements with a large parent group whose
    # mean and standard deviation is known.
    #
    # Input: a +Population+-descendant (typically a
    # +PopulationAPriori+) and a +Sample+ (typically a
    # +SampleArithmetic+)
    #
    # Output: a z-statistic, reflecting the probability that +sample+
    # was pulled from +population+.  To get a "proper" probability,
    # you need to run the z-statistic through +RStats.z_to_prob+ (or
    # use a table).
    #
    def self.zM(population, sample)
      n = sample.order
      m_s = sample.mean
      m_p = population.mean
      s_p = population.stddev

      z = Math.sqrt(n)*(m_p-m_s).abs/s_p
      return z
    end


    # The Student's t test (due to Gosset).  Langley 160-5
    #
    # This is a modified form of the zM test, in which you don't know
    # the standard deviation for the population.
    #
    # +sample+ must be of order 3 or more to get a .95 or greater
    # probability, so this raises SampleTooSmall if the input sample
    # is smaller than three.
    #
    # This function returns the t value, which is very similiar to the
    # z-statistic.  RStats does not currently implement a t-to-prob
    # function, so you'll need a table to interpret your results.
    #
    def self.students_t(population, sample)
      raise SampleTooSmall unless sample.order >= 3

      n = sample.order
      m_s = sample.mean
      s_s = sample.stddev
      m_p = population.mean

      t = Math.sqrt(n)*(m_p-m_s).abs/s_s
      
      return t
    end

    # Internal utility function:
    # Returns the rank value for both samples given. If +inv+ is
    # +true+, it returns the inverted rank value.
    #
    def self.sum_of_ranks_rank(sampleA, sampleB, inv=false)
      n_a = sampleA.order
      n_b = sampleB.order

      data_a = sampleA.data
      data_b = sampleB.data

      counts_a = Hash.new { |h,k| h[k] = 0 }
      data_a.each { |n| counts_a[n] += 1 }

      counts_b = Hash.new { |h,k| h[k] = 0 }
      data_b.each { |n| counts_b[n] += 1 }

      total_data = (data_a+data_b).sort
      total_data.reverse! if(inv)

      ranks = Hash.new { |h,k| h[k] = [] }
      total_data.each_with_index { |n,i| ranks[n].push(i+1) }

      ranks.each { |n,a| ranks[n] = RStats::mean_arith(a) }

      rank_a = counts_a.keys.inject(0.0) { |s,n| s += counts_a[n]*ranks[n] }
      rank_b = counts_b.keys.inject(0.0) { |s,n| s += counts_b[n]*ranks[n] }

      return rank_a, rank_b
    end

    # Internal utility function:
    # The core of ranks-based tests.
    def self.sum_of_ranks_kernel(sampleA, sampleB, skip_inv=false)
      n_a = sampleA.order
      n_b = sampleB.order

      rank_a, rank_b = sum_of_ranks_rank(sampleA, sampleB, false)

      r, n_r = (rank_a < rank_b ? [rank_a,n_a] : [rank_b,n_b])

      if(!skip_inv && n_a != n_b && n_a < 20 && n_b < 20)
	revrank_a, revrank_b = sum_of_ranks_rank(sampleA, sampleB, true)

	if(revrank_a < r)
	  r = revrank_a
	  n_r = n_a
	elsif(revrank_b < r)
	  r = revrank_b
	  n_r = n_b
	end
      end
      
      #puts "sum_of_ranks_kernel: #{r}/#{n_r} #{n_a}/#{rank_a} #{n_b}/#{rank_b}"

      return [r,n_r]
    end

    # Wilcoxon's Sum of Ranks test.  Langley pp166-78
    #
    # Purpose: to compare 2 unmatched random samples of measurements,
    # such as two samples taken from different sources.
    #
    # Input: two +Sample+s, which are run through sum_of_ranks_kernel,
    # then pushed through the appropriate formula.
    #
    # Output: a z-statistic, the probability that they are drawn from
    # the same population (ie: any differences observed are due to
    # chance, not to being from different groups).
    #
    def self.sum_of_ranks(sampleA, sampleB)
      r,n_r = sum_of_ranks_kernel(sampleA, sampleB)
      n_a = sampleA.order
      n_b = sampleB.order

      z = (n_r*(n_a+n_b+1)-2*r) / Math.sqrt(n_a*n_b*(n_a+n_b+1)/3)

      return z
    end



    # Wilcoxon's Signed Ranks Test. Langley pp179-189
    #
    # Purpose: To compare 2 random samples of measurements which are
    # matched.  The matching may be acheived by
    #   a) pairing members of 2 sample groups, carefully, or
    #   symmetrically, or by splitting
    #
    #   b) using a single sample group with different treatments,
    #   different methods of testing, different observers, or on
    #   different occasions.
    #
    # Input: a +SampleAMatched+, matched arithmetic sample.  The
    # sample must be of order 6 to get meaningful results.
    #
    # Output: a z-statistic, the probability that the matched samples
    # are drawn from the same population (ie: any differences observed
    # are due to chance, not to being from different groups).
    #
    def self.signed_ranks(msample)
      a = msample.a
      b = msample.b

      d = Array.new
      b.each_with_index { |n,i|
	cur_d = n - a[i]
	d.push(cur_d) unless (cur_d.abs < 0.000001)
      }

      n = 0.0+d.length

      raise SampleTooSmall if(n < 6)

      d_pos = Array.new
      d_neg = Array.new

      d.each { |d_i| (d_i > 0 ? d_pos : d_neg).push(d_i.abs) }

      s_pos = RStats::SampleArithmetic.new(d_pos)
      s_neg = RStats::SampleArithmetic.new(d_neg)

      r,n_r = RStats::Tests.sum_of_ranks_kernel(s_pos, s_neg, false)

      z = (n*(n+1)/2-2*r)/(Math.sqrt(n*(n+1)*(2*n+1)/6))

      #puts "SR: #{n} #{r} #{n_r} : #{z}"
      return z
    end


    # Wilcoxon's Stratified Test.  Langley pp190-198
    #
    # Purpose: To compare 2 independent stratified random samples of
    # measurements which have comparable strata, and the same number
    # of measurements in both samples.
    #
    # Input: +msamples+, a matched sample with equal lengths, etc.
    # Typically a +SampleAMatched+ matched arithmetic sample.
    #
    # Output: a z-statistic, the probability that the matched samples
    # are drawn from the same population (ie: any differences observed
    # are due to chance, not to being from different groups).
    #
    def self.stratified(msamples)
      k = msamples.length
      ns = msamples.map { |m| m.order }

      ranks = msamples.map { |m| 
	s_a = RStats::SampleArithmetic.new(m.a)
	s_b = RStats::SampleArithmetic.new(m.b)
	sum_of_ranks_rank(s_a, s_b)
      }

      rank_a = ranks.inject(0.0) { |s,r_i| s += r_i[0] }
      rank_b = ranks.inject(0.0) { |s,r_i| s += r_i[1] }

      r = rank_a < rank_b ? rank_a : rank_b

      z_num = ns.inject(-2*r) { |s,n| s += n*(2*n+1) }
      z_den = ns.inject(0.0) { |s,n| s += n*n*(2*n+1) }

      #puts "Strat: #{r} #{rank_a} #{rank_b} #{z_num}/sqrt(#{z_den}/3)"

      z = z_num/Math.sqrt(z_den/3)

      return z
    end

    # Spearman's Correlation Test.  Langley pp199-211
    #

    # Purpose: To test for correlation between two measurable
    # characteristics.  The measurable characteristics may exist
    # either
    #   a) in each individual of a sample group (e.g. weight and
    #   height), or
    #
    #   b) at the same time (e.g. the temperature in England and the
    #   wearing of overcoats in Germany)
    #
    # Input: a matched sample, +msample+ (typically of
    # +SampleAMatched+), where the two features to be correlated are
    # matched across the samples.  That is, if you were considering
    # height vs weight, the first part of the msample would be all the
    # weights, and the second part would be the corresponding heights.
    #
    # Output: a Spearman correlation coefficient, signed.  The sign of
    # the correlation coefficient gives the type of correlation
    # (positive/negative).  The degree of correlation can be gotten as
    # a z-statistic by taking the returned value through
    # +RStat.spearman_to_z+ (as a z-statistic).  The z-stat can then
    # be turned int a probability via +RStat.z_to_prob+.
    #
    def self.spearman_correlation(msample)
      n = msample.order

      raise SampleTooSmall unless n >= 5

      a = msample.a
      b = msample.b

      a_ranked = RStats::array_to_rank(a)
      b_ranked = RStats::array_to_rank(b)

      d2 = Array.new
      a_ranked.each_with_index { |a_i, i| 
	d2[i] = (a_i - b_ranked[i])*(a_i - b_ranked[i])
      }
      
      d2sum = d2.inject(0.0) { |s,d2_i| s += d2_i }

      occ_a = Hash.new { |h,k| h[k] = 0 }
      a_ranked.each { |a_i| occ_a[a_i] += 1 }

      occ_b = Hash.new { |h,k| h[k] = 0 }
      b_ranked.each { |b_i| occ_b[b_i] += 1 }
      
      ties = Hash.new { |h,k| h[k] = 0.0 }

      occ_a.values.each { |a_i| ties[0.0+a_i] += 1 }
      occ_b.values.each { |b_i| ties[0.0+b_i] += 1 }

      bigT = 0.0
      ties.each { |x, t_x|
	next if x == 1
	bigT += t_x*(x*x*x-x)/12
      }
      
      d2t = d2sum + bigT

      # puts "\nSpearman: #{d2t}: #{d2sum}+#{bigT} #{n}"

      co_factor = 6*d2t/(n*n*n-n)

      return co_factor
    end

    # Kruskal and Wallis' Test.  Langley pp 212-21
    #
    # Purpose: To compare 3 or more unmatched random samples of
    # measurements
    #
    # This test is subject to many details; you should really read the
    # introductory notes in Langley.  Or just not use it, as it's an
    # approximate test, and rather superfluous in the age of
    # computers.
    #
    # Input: +samples+, an array of 3 or more +Sample+s, typically all
    # +SampleArithmetic+
    #
    # Output: a chi^2 statistic reflecting the probability that the
    # variability between samples is due to chance (ie: they're not
    # from the same population).  This may be converted to a
    # probability with +RStat.pochisq+.
    #
    def self.kruskal_wallis(samples)
      ns = samples.map { |s| 0.0+s.order }

      bigN = ns.inject(0) { |s,n_i| s += n_i }

      data = samples.map { |s| s.data }
      all_data = data.inject([]) { |s,d_i| s = s + d_i }.sort
      all_ranks = RStats.array_to_rank(all_data)
      ranks = Hash.new
      all_data.each_with_index { |d_i, i| ranks[d_i] = all_ranks[i] }

      sample_ranks = data.map { |d| 
	d.inject(0.0) { |s,d_i| s += ranks[d_i] }
      }

      sr_sum = Array.new
      sample_ranks.each_with_index { |r_i,i| sr_sum[i] = r_i*r_i/ns[i] }
      ss = sr_sum.inject(0.0) { |s,r_i| s += r_i }

      #puts "\nKW: #{bigN} #{ns.inspect} #{sample_ranks.inspect} #{sr_sum.inspect} #{ss}"

      chi2 = 12.0/(bigN*bigN+bigN)*ss - 3*(bigN+1)
      return chi2
    end


    # Friedman's Test.  Langley pp222-9
    #
    # Purpose: To compare 3 or more random samples of measurements
    # which are matched.  The matching may be achieved by:
    #   a) carefully replicating members of each sample group, e.g. by
    #   dividing homogeneous substances or speciments into protions
    #   and allotting one portion to each sample group for different treatments, or
    #
    #   b) using 1 sample group for different treatments, different
    #   methods of testing, different observers, or on different
    #   occasions
    #
    # Input: +samples+, an +Array+ of co-indexed (matched) +Sample+s.
    #
    # Output: a chi^2 statistic reflecting the probability that the
    # variability between samples is due to chance (ie: they're not
    # from the same population).  This may be converted to a
    # probability with +RStat.pochisq+.
    #
    def self.friedman(samples)
      k = samples.length
      n = samples[0].order

      data = samples.map { |s| s.data }
      ranks = samples.map { 0.0 }

      n.times { |i|
	d_i = data.map { |d| d[i] }
	r = RStats::array_to_rank(d_i)
	r.each_with_index { |r_i, j| ranks[j] += r_i }
      }

      r = ranks.inject(0.0) { |s,r_i| s += r_i*r_i }
      chi2 = 12*r/(n*k*(k+1)) - 3*n*(k+1)

      return chi2
    end


    # Chi-Square contingency test.  Bass; Langley pp269-84
    #
    # Purpose: to check for "evidence of association between two
    # qualities (such as brains and beauty), when at least one of
    # these qualities is divided into three or more categories (such
    # as beautiful, average, and ugly)." [Langley, p269]
    #
    # Input: +data+ an +Array+ of +Array+s, representing the input
    # data matrix (row-major).
    #
    # Output: the chi^2 statistic reflecting the probability that the
    # qualities are not associated.
    #
    def self.chi_square_cont(data)
      row_sums = data.map { |row| row.inject(0.0) { |s,n| s += n } }

      n_rows = data[0].length

      col_sums = Array.new
      n_rows.times { |i|
	col_sums[i] = data.inject(0.0) { |s,r| s += r[i] }
      }

      bigN = row_sums.inject(0.0) { |s,r| s += r }

      row_probs = row_sums.map { |s| s/bigN }
      col_probs = col_sums.map { |s| s/bigN }

      chi2 = 0.0

      data.each_with_index { |row, i|
	row.each_with_index { |got, j|
	  f_e = row_probs[i]*col_probs[j]*bigN
	  #puts "chi2: #{i},#{j}: #{got} - #{f_e} (#{row_probs[i]}*#{col_probs[j]}"
	  chi2 += (got-f_e)*(got-f_e)/f_e
	}
      }
      return chi2      
    end



    # Chi-Square Goodness-of-fit test.  Bass; Langley pp269-84
    #
    # Purpose: much the same as zM above, to check if the variations
    # between two samples of measurements are statistically
    # significant.
    #
    # Input: +before+, one sample (as an array of measurements),
    # +after+, the other sample (also as array of measurements).
    #
    # Output: the chi^2 statistic reflecting the probability that the
    # differences in measurements is due to chance.  This may be
    # converted to a proper probability with +RStats.pochisq+.
    #
    def self.chi_square_gof(before, after)
      # Normalize to expectations
      after_sum = after.inject(0.0) { |s,n| s += n }
      before_sum = before.inject(0.0) { |s,n| s += n }

      norm = after_sum/before_sum

      expected = before.map { |n| norm * n }

      # chi2 = sum (fa-fe)^2/fe

      s = 0.0
      after.each_with_index { |a_i, i| 
	s += (a_i-expected[i])*(a_i-expected[i])/expected[i]
      }

      return s
    end
  end # module Tests
end # module RStats
