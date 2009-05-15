#!/usr/bin/ruby -Ilib/ -w

# Run these from rstats/ as ruby ./test/tc_rstats.rb

require 'test/unit'
require 'rstats'

ARITH_TESTS = [
  [ 0.0, 0.0, [ 0.0, 0.0, 0.0 ] ],
  [ 0.0, 1.0, [ -1.0, 0.0, 1.0 ] ],
  [ 1.0, 0.0, [ 1.0, 1.0, 1.0 ] ],
  [ 3.0, 1.5811, [ 1, 2, 3, 4, 5 ] ],
  [ -1.0, 0.0, [ -1.0, -1.0, -1.0] ],
  [ -1.0, 1.0, [ -2.0, -1.0, 0.0] ],
]

GEO_TESTS = [
]


class TestStats < Test::Unit::TestCase
  def assert_within(a,b,d)
    puts "\n#{a} not within #{b}+/-#{d}" if((a-b).abs > d)
    assert( (a-b).abs <= d)
  end

  def test_mean_arith_assert()
    [nil, 42, []].each { |d|
      assert_raises(RStats::BogusData) {
	RStats::mean_arith(d)
      }
    }
  end

  def test_stddev_arith_asserts()
    [nil, 42, [] ].each { |d|
      assert_raises(RStats::BogusData) {
	RStats::stddev_arith(d, 1.0)
      }
    }
  end

  def test_mean_stddev_arith()
    ARITH_TESTS.each { |t|
      m = t[0]
      s = t[1]
      d = t[2]

      # puts "Check #{m}/#{s}/#{d.inspect} (#{d.length})"

      assert(RStats::mean_arith(d) == m)
      ss = RStats::stddev_arith(d,m)
      assert_within(ss, s, 0.0001)
    }
  end

  def test_sample_arith()
    ARITH_TESTS.each { |t|
      m = t[0]
      s = t[1]
      d = t[2]

      sample = RStats::SampleArithmetic.new(d)
      assert_within(sample.mean, m, 0.0001)
      assert_within(sample.stddev, s, 0.0001)
    }
  end

  def test_mean_geo_assert()
    [nil, 42, []].each { |d|
      assert_raises(RStats::BogusData) {
	RStats::mean_geo(d)
      }
    }
  end

  def test_stddev_geo_asserts()
    [nil, 42, [], [1]].each { |d|
      assert_raises(RStats::BogusData) {
	RStats::stddev_geo(d, 1.0)
      }
    }
  end

  def test_mean_stddev_geo()
    GEO_TESTS.each { |t|
      m = t[0]
      s = t[1]
      d = t[2]

      # puts "Check #{m}/#{s}/#{d.inspect} (#{d.length})"

      assert(RStats::mean_geo(d) == m)
      ss = RStats::stddev_geo(d,m)
      assert_within(ss, s, 0.0001) if s
      puts "got geo stddev #{ss}" if !s
    }
  end

  def test_sample_geo()
    GEO_TESTS.each { |t|
      m = t[0]
      s = t[1]
      d = t[2]

      sample = RStats::SampleGeometric.new(d)
      assert_within(sample.mean, m, 0.0001)
      assert_within(sample.stddev, s, 0.0001)
    }
  end


  def test_test_zm()
    tests = [
      #  m    n   M    S  z      P
      [ 1072, 1, 1060, 3, 4.00, 0.00013 ],
      [ 73,  40, 70,   5, 3.79, 0.00029 ],
      [ 0.637, 14, 0.744, 0.262, 1.53, 0.12412 ],
    ]

    tests.each { |t|
      m,n,p_m,p_s,z_given,prob = t

      pop = RStats::PopulationAPriori.new(p_m, p_s)

      data = [m] * n
      samp = RStats::SampleArithmetic.new(data)

      z = RStats::Tests.zM(pop, samp)
      assert_within(z_given, z, 0.01)
      p = RStats::z_to_prob(z)
      assert_within(prob, p, 0.00001)
    }
  end


  def test_test_zm_copper()
    # Langly, p158&9
    # Copper melts at average of 1080degC, with s=5degC
    cuPop = RStats::PopulationAPriori.new(1080, 5)

    # Given a sample that melts at 1072, is it likely to be copper?
    sample_single = RStats::SampleArithmetic.new([1072])
    z_single = RStats::Tests.zM(cuPop, sample_single)
    p_single = RStats::z_to_prob(z_single)
    #puts "p_single: #{p_single}"

    # This should be ~0.11, so there's about a 1 in 9 chance of it
    # being copper
    assert_within(p_single, 0.1109, 0.0001)

    # Now, try it three more times, getting 1071,2,3:
    sample_four = RStats::SampleArithmetic.new([1072,1071,1072,1073])
    z_four = RStats::Tests.zM(cuPop, sample_four)
    p_four = RStats::z_to_prob(z_four)
    #puts "p_four: #{p_four}"

    # This probability is less than 0.2%, so it's almost certain
    # that the four samples are not from the copper population.
    assert_within(p_four, 0.0023, 0.0001)
  end

  def test_zm_haircuts()
    london = RStats::PopulationAPriori.new(0.954, 0.600)
    paris = RStats::SampleArithmetic.new([0.764]*40)

    z_same_pop = RStats::Tests.zM(london, paris)
    p_same_pop = RStats::z_to_prob(z_same_pop)
    #puts "There is a #{p_same_pop} chance that the difference is not significant." 
    assert_within(p_same_pop, 0.0536, 0.0001)
  end

  def test_studentst_press()
    pop = RStats::PopulationAPriori.new(45, 0)
    sample = RStats::SampleArithmetic.new([46,47,48])

    t = RStats::Tests.students_t(pop, sample)
    assert_within(t, 3.46, 0.01)


    sample = RStats::SampleArithmetic.new([46,47,48,47,47])
    t = RStats::Tests.students_t(pop, sample)
    assert_within(t, 6.32, 0.01)
  end

  def test_studentst_phone()
    pop = RStats::PopulationAPriori.new(48, 0)
    sample = RStats::SampleArithmetic.new([56, 51, 63, 60])

    t = RStats::Tests.students_t(pop, sample)
    assert_within(t, 3.65, 0.01)
  end

  def test_sum_of_ranks_op()
    sa = RStats::SampleArithmetic.new([16, 20, 25, 19, 22, 15, 22, 19])
    sb = RStats::SampleArithmetic.new([18, 19, 15, 16, 21, 17, 17, 14])

    z = RStats::Tests.sum_of_ranks(sa,sb)
    p = RStats::z_to_prob(z)
    assert_within(p, 0.096, 0.001)
  end

  def test_sum_of_ranks_ads()
    s_who = RStats::SampleArithmetic.new([341, 326, 3260, 305, 326])
    s_why = RStats::SampleArithmetic.new([352, 382, 347])

    z = RStats::Tests.sum_of_ranks(s_who, s_why)
    p = RStats::z_to_prob(z)
    assert_within(p, 0.162, 0.001)
  end

  def test_sum_of_ranks_mass()
    s_fav = RStats::SampleArithmetic.new([3.05, 3.01, 3.20, 3.16, 3.11, 3.09])
    s_new = RStats::SampleArithmetic.new([3.18, 3.23, 3.19, 3.28, 3.08, 3.18])

    z = RStats::Tests.sum_of_ranks(s_fav, s_new)
    p = RStats::z_to_prob(z)
    assert_within(p, 0.084, 0.001)
  end

  def test_sum_of_ranks_routes()
    s_old = RStats::SampleArithmetic.new([3204, 2967, 3053, 3267, 3370, 3492, 3105, 3330])
    s_new = RStats::SampleArithmetic.new([3568, 3299, 3618, 3494])

    z = RStats::Tests.sum_of_ranks(s_old, s_new)
    p = RStats::z_to_prob(z)
    assert_within(p, 0.0344, 0.0001)
  end

  def test_signed_ranks_nockout
    t_pheno = [ 7.5, 7, 7   , 5.75, 4.25, 9.25, 8  , 7.25, 8.5, 7.75]
    t_nock  = [ 8  , 6, 6.75, 5,    4.5 , 8   , 7.5, 6.25, 8  , 7.75]

    ms = RStats::SampleAMatched.new([t_pheno, t_nock])
    z = RStats::Tests.signed_ranks(ms)
    p = RStats::z_to_prob(z)

    assert_within(p, 0.052, 0.001)
  end

  def test_signed_ranks_hair()
    t_pre = [105, 105, 93, 120, 111, 80, 91]
    t_bleach = [ 97, 95, 93, 117, 108, 85, 86]

    ms = RStats::SampleAMatched.new([t_pre, t_bleach])
    z = RStats::Tests.signed_ranks(ms)
    p = RStats::z_to_prob(z)

    assert_within(p, 0.136, 0.001)
  end

  def test_signed_ranks_hair2()
    t_pre = [105, 105, 93.2, 120.1, 111.4, 80.1, 91.3]
    t_bleach = [ 97, 95, 93.0, 117.1, 108.3, 84.7, 86.0]

    ms = RStats::SampleAMatched.new([t_pre, t_bleach])
    z = RStats::Tests.signed_ranks(ms)
    p = RStats::z_to_prob(z)

    assert_within(p, 0.0956, 0.001)
  end

  def test_signed_ranks_car()
    d_gas = [17.1, 29.5, 23.8, 37.3, 19.6, 24.2, 30.0, 20.9]
    d_lub = [14.2, 30.3, 21.5, 36.3, 19.6, 24.5, 26.7, 20.6]

    ms = RStats::SampleAMatched.new([d_gas, d_lub])
    z = RStats::Tests.signed_ranks(ms)
    p = RStats::z_to_prob(z)
    assert_within(p, 0.125, 0.001)
  end

  def test_signed_ranks_darwin()
    d_cross = [23.5, 12, 21,22, 19.125, 21.5, 22.125, 20.375, 18.25, 21.675, 23.25, 21, 22.125, 23, 12]
    d_self  = [17.375, 20.375, 20, 20, 18.375, 18.675, 18.675, 15.25, 16.5, 18, 16.25, 18, 12.75, 15.5, 18]

    ms = RStats::SampleAMatched.new([d_cross, d_self])
    z = RStats::Tests.signed_ranks(ms)
    p = RStats::z_to_prob(z)
    assert_within(p, 0.00509, 0.0001)
  end

  def test_sum_of_ranks_sleep()
    # This is just here to answer Q52, Langley p189
    t_pheno = [ 7.5, 7, 7   , 5.75, 4.25, 9.25, 8  , 7.25, 8.5, 7.75]
    t_nock  = [ 8  , 6, 6.75, 5,    4.5 , 8   , 7.5, 6.25, 8  , 7.75]

    s_pheno = RStats::SampleArithmetic.new(t_pheno)
    s_nock =  RStats::SampleArithmetic.new(t_nock)

    z = RStats::Tests.sum_of_ranks(s_pheno, s_nock)
    p = RStats::z_to_prob(z)

    #puts "\nSum-of-ranks pheno: #{z} #{p}"

    ms = RStats::SampleAMatched.new([t_pheno, t_nock])
    z = RStats::Tests.signed_ranks(ms)
    p = RStats::z_to_prob(z)

    #puts "Signed ranks phen: #{z} #{p}"
  end

  def test_test_stratified_ddt()
    data = [
      [ [ 18, 26, 30, 50 ], [ 34, 42, 53, 63 ] ],
      [ [ 33, 42, 44, 44 ], [ 60, 62, 66, 80 ] ],
      [ [ 44, 50, 56, 64 ], [ 74, 77, 84, 92 ] ],
    ]

    msamples = data.map { |pair| RStats::SampleAMatched.new(pair) }

    z = RStats::Tests.stratified(msamples)
    p = RStats::z_to_prob(z)
    assert_within(p, 0.0005, 0.0001)
  end

  def test_test_stratified_acne()
    data = [
      [ [2, 3], [2, 4] ],
      [ [3, 5, 6, 10], [4, 6, 7, 9] ],
      [ [6, 8, 11], [9, 14, 14] ],
      [ [8, 10, 11], [12, 14, 15] ],
    ]
    msamples = data.map { |pair| RStats::SampleAMatched.new(pair) }

    z = RStats::Tests.stratified(msamples)
    assert_within(z, 2.03, 0.01)
    p = RStats::z_to_prob(z)
    assert_within(p, 0.050, 0.001)
  end

  def test_test_stratified_houdini()
    data = [
      [ [ 10.9, 11.3, 10.2 ], [ 10.3, 10.6, 12.2 ] ],
      [ [ 13.8, 15.1, 14.3 ], [ 16.3, 15.2, 15.8 ] ],
    ]

    msamples = data.map { |pair| RStats::SampleAMatched.new(pair) }

    z = RStats::Tests.stratified(msamples)
    p = RStats::z_to_prob(z)
    assert_within(p, 0.121, 0.001)
  end

  def test_test_stratified_powder()
    data = [
      [ [0,1], [3,3] ],
      [ [2,5], [4,6] ],
      [ [6,4], [5,7] ],
      [ [10,8], [7,11] ],
      [ [4,1], [3,2] ],
      [ [1,2], [1,4] ],
      [ [6,5], [9,9] ],
      [ [0,2], [2,3] ],
      [ [4,7], [3,5] ]
    ]

    msamples = data.map { |pair| RStats::SampleAMatched.new(pair) }

    z = RStats::Tests.stratified(msamples)
    assert_within(z, 1.81, 0.01)
    p = RStats::z_to_prob(z)

    #puts "Got powder: #{z} #{p}"
    assert_within(p, 0.078, 0.001)
  end

  def test_test_spearman_auradio()
    radio = [ 171, 178, 251, 160, 155 ]
    tele  = [  74, 224, 300, 404, 323 ]

    msample = RStats::SampleAMatched.new([radio,tele])

    co_factor = RStats::Tests.spearman_correlation(msample)
    z = RStats.spearman_to_z(msample.order, co_factor)
    p = RStats.z_to_prob(z)

    # puts "\nSpearman TV: got #{co_factor} -> #{z} -> #{p}"
    assert_within(p, 0.241, 0.001)
  end

  def test_test_spearman_beauty()
    ages   = [ 17, 16, 18, 20, 18, 18, 20, 23 ]
    places = [  1,  2,  2,  4,  5,  6,  7,  8 ]

    msample = RStats::SampleAMatched.new([ages, places])

    co_factor = RStats::Tests.spearman_correlation(msample)
    z = RStats.spearman_to_z(msample.order, co_factor)
    p = RStats.z_to_prob(z)

    #puts "\nSpearman Beauty: got #{co_factor} -> #{z} -> #{p}"
    assert_within(p, 0.052, 0.001)
  end

  def test_test_spearman_spelling()
    # This is a test correlating foot length (inches) to spelling
    # scores... Langley pp208-10

    feet = [ 6.5, 6.75, 6.75, 7, 7.25, 7.5, 7.5, 7.5, 7.5, 7.75, 8, 8,
    8, 8, 8.25, 8.5, 8.75, 9, 9.25, 9.25, 9.5, 9.5, 9.5, 9.5, 9.5,
    9.75, 9.75, 10, 10, 10.25, 10.25, 10.5, 10.5, 10.5, 10.5, 10.75,
    10.75 ]

    scores = [ 16, 28, 46, 14, 41, 10, 56, 43, 15, 21, 50, 28, 57, 65,
    42, 36, 71, 47, 47, 66, 71, 71, 61, 55, 86, 62, 78, 71, 59, 60,
    63, 68, 98, 88, 91, 98, 78]

    msample = RStats::SampleAMatched.new([feet, scores])

    co_factor = RStats::Tests.spearman_correlation(msample)
    z = RStats.spearman_to_z(msample.order, co_factor)
    assert_within(z, 4.94, 0.01)
    p = RStats.z_to_prob(z)
    assert_within(p, 0.00, 0.001)
  end

  def test_test_spearman_chem()
    temp = [ 15, 18, 21, 24, 27, 30, 33 ]
    cyield = [ 66, 69, 69, 70, 64, 73, 75 ]
    msample = RStats::SampleAMatched.new([temp, cyield])
    
    co_factor = RStats::Tests.spearman_correlation(msample)
    z = RStats.spearman_to_z(msample.order, co_factor)
    p = RStats.z_to_prob(z)
    assert_within(p, 0.123, 0.001)
    #puts "\nSpearman chemical yield: #{co_factor} -> #{z} -> #{p}"
  end
  
  def test_test_spearman_aupost()
    letters = [ 1486, 1448, 1473, 1570, 1619, 1705 ]
    offices = [ 8315, 8315, 8261, 8244, 8234, 8222 ]
    msample = RStats::SampleAMatched.new([letters, offices])
    
    co_factor = RStats::Tests.spearman_correlation(msample)
    z = RStats.spearman_to_z(msample.order, co_factor)
    assert_within(z, 1.98, 0.01)
    p = RStats.z_to_prob(z)
    assert_within(p, 0.056, 0.001)
    #puts "\nSpearman AU post: #{co_factor} -> #{z} -> #{p}"
  end
 
  def test_test_krusal_wallis_bottlecaps()
    data = [ 
      [340, 345, 330, 342, 338],
      [339, 333, 344],
      [347, 343, 349, 355]
    ]

    samples = data.map { |d| RStats::SampleArithmetic.new(d) }

    chi2 = RStats::Tests.kruskal_wallis(samples)
    assert_within(chi2, 5.66, 0.01)
    p = RStats.pochisq(chi2, 1)

    #puts "\nKW bottlecaps: got #{chi2}, #{p}"
  end


  def test_test_krusal_wallis_singers()
    data = [
      [36, 22, 19, 16],
      [39, 14, 20, 18],
      [21, 32, 28, 22],
    ]

    samples = data.map { |d| RStats::SampleArithmetic.new(d) }

    chi2 = RStats::Tests.kruskal_wallis(samples)
    assert_within(chi2, 1.3, 0.1)
    p = RStats.pochisq(chi2, 3)

    #puts "\nKW singers: got #{chi2}, #{p}"
    
  end

  def test_test_friedman_golf_player()
    data = [
      [80, 80, 85, 90, 85, 81],
      [77, 81, 82, 86, 80, 82],
      [81, 83, 84, 85, 86, 82],
      [82, 85, 87, 87, 81, 79],
    ]

    samples = data.map { |d| RStats::SampleArithmetic.new(d) }

    chi2 = RStats::Tests.friedman(samples)
    assert_within(chi2, 3.15, 0.01)
    # puts "\nfm golf play: #{chi2}"
    p = RStats.pochisq(chi2, 4)

    #puts "\nFriedman golf players got #{chi2}, #{p}"

  end


  def test_test_friedman_golf_times()
    data = [
      [80, 77, 81, 82],
      [80, 81, 83, 85],
      [85, 82, 84, 87],
      [90, 86, 85, 87],
      [85, 80, 86, 81],
      [81, 82, 82, 79],
    ]

    samples = data.map { |d| RStats::SampleArithmetic.new(d) }

    chi2 = RStats::Tests.friedman(samples)
    assert_within(chi2, 11.96, 0.01)
    #puts "\nfm golf time: #{chi2}"
    p = RStats.pochisq(chi2, 4)

    #puts "\nFriedman golf times: got #{chi2}, #{p}"
  end

  def test_test_friedman_rolls()
    data = [
      [ 3, 2, 2, 3, 3  , 2, 3, 1],
      [ 2, 3, 3, 1, 1.5, 3, 1, 3],
      [ 1, 1, 1, 2, 1.5, 1, 2, 2]
    ]

    samples = data.map { |d| RStats::SampleArithmetic.new(d) }

    chi2 = RStats::Tests.friedman(samples)
    assert_within(chi2, 3.94, 0.01)
    #puts "\nfm rolls: #{chi2}"
  end

  def test_test_chi2gof()
    before = [ 5,  9 , 65, 10, 6, 5]
    after  = [22, 45, 198, 30, 9, 1]

    chi2 = RStats::Tests.chi_square_gof(before, after)
    assert_within(chi2, 32.26, 0.01)
  end


  def test_test_chi2_cont()
    data = [ [79,48],[1091,1492] ]
    chi2 = RStats::Tests.chi_square_cont(data)
    assert_within(chi2, 19.67, 0.01)
  end

  def test_test_chi2_gof_accidents()
    accidents = [ 7, 7, 1 ]
    expected = [ 1, 1, 1 ]
    chi2 = RStats::Tests.chi_square_gof(expected, accidents)
    assert_within(chi2, 4.80, 0.01)
  end

  def test_test_chi2_gof_aupop()
    expected = [30, 20, 21, 16, 10, 3]
    patients = [19, 32, 22, 19,  6, 2]
    chi2 = RStats::Tests.chi_square_gof(expected, patients)
    assert_within(chi2, 13.77, 0.01)
  end

  def test_test_chi2_cont_cholera()
    data = [ [3, 276], [66, 473] ]
    chi2 = RStats::Tests.chi_square_cont(data)
    #puts "\nchi2: #{chi2}"
  end

  def test_test_chi2_cont_sunlight_eyes()
    data = [ 
      [19, 27, 4],
      [ 7,  8, 5],
      [ 1, 13, 16]
    ]
    chi2 = RStats::Tests.chi_square_cont(data)
    # puts "\nchi2: #{chi2}"
    assert_within(chi2, 25.13, 0.01)
  end

  def test_test_chi2_cont_exam()
    data = [
      [ 10, 45,  5 ],
      [  4, 35, 11 ],
    ]
    chi2 = RStats::Tests.chi_square_cont(data)
    assert_within(chi2, 5.20, 0.01)
  end

  def test_test_chi2_gof_ski()
    expected = [ 1, 1, 1]
    accidents = [6, 11, 16]

    chi2 = RStats::Tests.chi_square_gof(expected, accidents)
    assert_within(chi2, 4.55, 0.01)
  end  

  def test_test_chi2_gof_ps19_3()
    expected = [ 20, 20, 20, 20, 20, 20 ]
    rolls = [ 18, 23, 16, 21, 18, 24 ]

    chi2 = RStats::Tests.chi_square_gof(expected, rolls)
    assert_within(chi2, 2.5, 0.1)
  end

  def test_test_chi2_cont_ps19_4()
    data = [
      [ 50, 87,  5,  8 ],
      [ 40, 69, 60, 11 ],
      [ 15, 13, 42,  5 ],
      [  5, 27, 17,  1 ],
      [ 15,  4,  1, 25 ],
    ]

    chi2 = RStats::Tests.chi_square_cont(data)
    assert_within(chi2, 223, 1) # XXX 217 in the book
  end

  def test_test_chi2_cont_ps19_5()
    data = [
      [ 25,  5 ],
      [ 10,  5 ],
      [  5, 10 ],
    ]

    chi2 = RStats::Tests.chi_square_cont(data)
    assert_within(chi2, 11.25, 0.01)
  end

  def test_test_chi2_cont_ps19_6()
    data = [
      [ 10,  6, 13, 13 ],
      [ 10, 12, 19, 21 ],
      [ 13, 10, 13, 18 ]
    ]

    chi2 = RStats::Tests.chi_square_cont(data)
    assert_within(chi2, 2.13, 0.01)
  end
  
  def test_test_chi2_cont_ps19_7()
    data = [
      [ 10, 10, 10, 20 ],
      [ 10, 40, 50, 50 ],
      [ 13, 25, 60, 52 ],
      [ 17, 25, 30, 78 ],
    ]

    chi2 = RStats::Tests.chi_square_cont(data)
    assert_within(chi2, 32.56, 0.01)
  end


end

