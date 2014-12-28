# A bunch of statistics stuff in ruby.
#
# Copyright (c) 2006 Josh Myer <josh@joshisanerd.com>

# Mostly just a programmatic form of notes from Russell Langley's
# "Practical Statistics Simply Explained", Dover.  ISBN 0486227294

#  The book is highly recommended as a primer, though it is a little
#  light on the theoretical underpinnings. 
# A bunch of statistics stuff in ruby.
#
# Copyright (c) 2006 Josh Myer <josh@joshisanerd.com>

# Mostly just a programmatic form of notes from Russell Langley's
# "Practical Statistics Simply Explained", Dover.  ISBN 0486227294

#  The book is highly recommended as a primer, though it is a little
#  light on the theoretical underpinnings. 

require 'rstats'


module RStats

  # A sample of data distributed such that arithmetic methods (mean,
  # stddev, etc) make sense.
  class SampleArithmetic
    attr_reader(:data, :mean, :stddev)
    def initialize(points)
      @data = points
      @mean = RStats::mean_arith(@data)
      @stddev = RStats::stddev_arith(@data, @mean)
    end

    def order()
      return @data.length
    end
  end

  class SampleGeometric < SampleArithmetic
    def initialize(points)
      @data = points
      @mean = RStats::mean_geo(@data)
      @stddev = RStats::stddev_geo(@data, @mean)
    end
  end

  # A matched sample, distributed arithmetically.
  #
  class SampleAMatched < SampleArithmetic
    # +mpoints+ is an array of two arrays, one for each set of data
    def initialize(mpoints)
      raise BogusData if(!mpoints || !mpoints.respond_to?("length") || !mpoints.length == 2)

      [0,1].each { |i|
	raise BogusData if(!mpoints[i] || !mpoints[i].respond_to?("length"))
      }

      raise MismatchedSample unless mpoints[0].length == mpoints[1].length

      @data = mpoints
    end

    # Get the first (zero-th) set of data points
    def a()
      return @data[0]
    end

    # Get the second (one-th) set of data points
    def b()
      return @data[1]
    end

    def order()
      return @data[0].length
    end

  end # class SampleAMatched
end # module RStats
