require File.expand_path(__FILE__).sub(%r(/www/.*), '/test/test_helper.rb')
require File.expand_path(__FILE__).sub(%r(.*/test/), '').sub(/test_(.*)\.rb/,'\1')

require 'rbbt-util'

class TestHelpers < Test::Unit::TestCase
  def test_data_servers
    assert Array === Covid19.data_servers
  end

  def test_sample_info
    iii Covid19.sample_info
  end
end

