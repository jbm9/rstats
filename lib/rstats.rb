# A bunch of statistics stuff in ruby.
#
# Copyright (c) 2006 Josh Myer <josh@joshisanerd.com>
#
# Mostly just a programmatic form of notes from Russell Langley's
# "Practical Statistics Simply Explained", Dover.  ISBN 0486227294

#  The book is highly recommended as a primer, though it is a little
#  light on the theoretical underpinnings. 

# Introduce our namespace, just in case
module RStats
end # module RStats

require 'rstats/base'   # The Stats.foo() methods, Population
require 'rstats/sample' # The Sample classes
require 'rstats/tests'  # The Stats::Tests methods

