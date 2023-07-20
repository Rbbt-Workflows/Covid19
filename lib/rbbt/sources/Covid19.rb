require 'rbbt-util'
require 'rbbt/resource'

module Covid19
  extend Resource
  self.subdir = 'share/databases/Covid19'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  Covid19.claim Covid19.GSE145926, :proc do |directory|
    url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE145926\&format=file"
    Misc.untar(Open.open(url), directory)
    nil
  end


  Covid19.claim Covid19.models, :proc do |directory|
    Open.mkdir directory
    Rbbt.data.find(:lib).glob("*").each do |file|
      FileUtils.cp_r file, directory
    end
  end
end
