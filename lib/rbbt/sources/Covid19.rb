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

  Covid19.claim Covid19.models.epithelial_cell_2, :proc do |directory|
    Open.mkdir directory
    Rbbt.data.glob("epithelial_cell_2.*").each do |file|
      Open.cp file, File.join(File.dirname(directory), File.basename(file))
    end
    nil
  end
end
