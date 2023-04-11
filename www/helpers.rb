require 'rbbt/util/ssh'

module Covid19

  SERVERS = %w(turbo mn1)

  def self.persist(name, type, &block)
    Persist.persist(name, type, :persist_dir => 'var/cache/persist/Covid19', &block) 
  end

  def self.line(server)
    SSHLine.new(server)
  end

  def self.data_servers
    persist("data_servers", :array) do
      SERVERS.select{|s| SSHLine.workflow(s, "Covid19", "Covid19.glob('*').any?") }
    end
  end
end

