# A bunch of statistics stuff in ruby.
#
# Copyright (c) 2006 Josh Myer <josh@joshisanerd.com>

# Mostly just a programmatic form of notes from Russell Langley's
# "Practical Statistics Simply Explained", Dover.  ISBN 0486227294

#  The book is highly recommended as a primer, though it is a little
#  light on the theoretical underpinnings. 


require 'rstats'
BIGX = 20.0

module RStats
  class BogusData < Exception; end
  class BogusPopulation < Exception; end
  class SampleTooSmall < Exception; end
  class MismatchedSample < Exception; end

  INV_ROOT_2PI = 1.0 / Math.sqrt(2*Math::PI)

  # Find the arithmetic mean of a given array, +data+
  #
  def self.mean_arith(data)
    raise BogusData unless data && 
      data.respond_to?("length") && 
      data.length != 0

    sum = data.inject(0.0) { |sum,d| sum += d }

    return sum/data.length
  end

  # Find the standard deviation of a data set, given as an array
  # (+data+) and a mean (+mean+)
  #
  def self.stddev_arith(data, mean)
    raise BogusData unless data && 
      data.respond_to?("length") && 
      data.length > 0

    return 0 unless data.length > 1

    ss = data.inject(0.0) { |s,d| s += (d-mean)*(d-mean) }

    return Math.sqrt(ss/(data.length-1))
  end


  def self.mean_geo(data)
    raise BogusData unless data && 
      data.respond_to?("length") && 
      data.length != 0
    
    sum = data.inject(0.0) { |sum,d| sum += Math.log(d) }

    return Math.exp(sum/data.length)
  end

  def self.stddev_geo(data, mean)
    raise BogusData unless data && 
      data.respond_to?("length") && 
      data.length >= 2

    ss = data.inject(0.0) { |s,d| s += (Math.log(d)-mean)*(Math.log(d)-mean) }

    return Math.exp(Math.sqrt(ss/(data.length-1)))
  end

  # The template for a Population, not currently used for much more
  # than a placeholder for a mean and standard deviation.
  #
  # All methods of +Population+ raise +BogusPopulation+ unless you've
  # written a child class that overrides the method (think of this as
  # a pure virtual class)
  #
  class Population
    def initialize()
      @mean = nil
      @stddev = nil
    end

    def mean()
      raise BogusPopulation
    end

    def stddev()
      raise BogusPopulation
    end
  end

  # An a priori distribution, for which we know the mean and standard
  # deviation
  class PopulationAPriori < Population
    attr_reader(:mean, :stddev)

    def initialize(mean, stddev)
      @mean = mean
      @stddev = stddev
    end
  end # class PopulationAPriori


  # Convert a z-statistic to a probability
  def self.z_to_prob(z)
    return 1/Math.sqrt(2*Math::PI)*Math.exp(-z*z/2)
  end


  # Map an array to the rank-value of each item.  That is, the array
  # [ 8, 22, 18.5 ]
  # returns
  # [ 1, 3, 2 ]
  #
  # This is commonly used as a form of normalization
  #
  def self.array_to_rank(a)
    rank = Hash.new { |h,k| h[k] = Array.new }
    a.sort.each_with_index { |n,i| rank[n].push(i+1) }
    rank.each { |n,a_n| rank[n] = RStats::mean_arith(a_n) }
    return a.map { |n| rank[n] }
  end

  # Convert a Spearman coefficient to a z-statistic
  def self.spearman_to_z(n, co_factor)
    return Math.sqrt(n-1)*(1-co_factor).abs
  end

  # An exponential function, as found in the public-domain pochisq
  def self.ex(x) (x < -BIGX ? 0.0 : Math.exp(x)); end

  # probability-of-chi2, after the public-domain implementation in C
  # +x+ is the chi-squared value, +df+ is the degrees of freedom
  def self.pochisq(x, df)
    return 1.0 if(x <= 0.0 || df < 1)

    a = y = s = 0.0
    e = c = z = 0.0

    even = false

    a = 0.5 * x
    even = ((2*(df/2))==df)

    y = ex(-a) if df > 1

    s = (even ? y : 2.0*RStats::z_to_prob(-Math.sqrt(x)))

    if(df > 2)
      x = 0.5 * (df - 1.0)
      z = (even ? 1.0 : 0.5)
      if(a > BIGX)
	e = even ? 0.0 : Math.log(Math.sqrt(MATH::PI))
	c = Math.log(a)
	while(z <= x)
	  e = Math.log(z) + e
	  s += math.exp(c*z-a-e)
	  z += 1.0
	end
	return s;
      else
	e = (even ? 1.0 : 1/Math.sqrt(a*Math::PI))
	c = 0.0
	while(z <= x)
	  e = e * (a/z)
	  c = c + e
	  z += 1.0
	end
	return (c * y + s)
      end
    else
      return s
    end
  end

end # module RStats
