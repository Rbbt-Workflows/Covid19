Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/Covid19'

Workflow.require_workflow "PerMedCoE"

module Covid19
  extend Workflow


  desc <<~EOF
  Known sample info in GSE145926 dataset
  EOF
  task :sample_info => :tsv do
    file = Rbbt.data.metadata
    tsv = file.tsv :header_hash => "", :type => :list
    tsv.delete_if{|k,v| k == 'C146' }
    tsv.key_field = "id"
    tsv.fields = ["group", "file"]
    tsv.process "file" do |file|
      Covid19.GSE145926[File.basename(file)]
    end
    tsv
  end

  desc <<~EOF
  Process single cell data for a patient
  EOF
  dep_task :single_cell_processing, PerMedCoE, :single_cell_processing_BB, 
    :p_id => :jobname,
    :p_file => nil do |jobname,options|

      options[:p_id] = jobname

      if options[:p_file].nil?
        file = Covid19.GSE145926.produce.glob("*_#{jobname}_*.h5").first 
        #file = Rbbt.identify(file)
        raise ParameterException, "File not found for sample id #{jobname} (jobname)" if file.nil?
        options[:p_file] = Path.setup(file)
      end

      {:inputs => options}
    end


  desc <<~EOF
  Use the single cell processing results in the MaBoSS knockouts to
  personalize a patient
  EOF
  dep PerMedCoE, :MaBoSS_BB, :jobname => "Default",
    :positional => 'default', :model => 'epithelial_cell_2', :data_folder => Covid19.models,
    :model_folder => nil, :genes_druggable => nil, :genes_target => nil
  dep :single_cell_processing
  dep_task :personalize_single_cell_patient, PerMedCoE, :personalize_patient_BB, 
    :model_prefix => Covid19.models.epithelial_cell_2, :t => "T", 
    :positional => 'default', :expression => nil, :cnv => nil, :mutation => nil, :cell_type => nil,:model_bnd => nil, :model_cfg => nil,
    :ko => :placeholder, :norm_data => :placeholder, :cells => :placeholder  do |jobname,options,dependencies|

      maboss, single_cell = dependencies.flatten
      options[:norm_data] = single_cell.file('output').norm_data
      options[:cells] = single_cell.file('output').cells_metadata
      options[:ko] = maboss.file('output').ko_file
      {:inputs => options}
    end


  desc <<~EOF
  Run PhysiBoSS using the personalized models
  EOF
  dep :personalize_single_cell_patient
  dep PerMedCoE, :PhysiBoSS_BB,
    :bnd_file => :placeholder, :cfg_file => :placeholder, :sample => :placeholder, :prefix => :placeholder  do |jobname,options,dependencies|

      patient = dependencies.flatten.first
      maboss = patient.step(:MaBoSS_BB)
      options[:sample] = jobname
      options[:prefix] = jobname

      personalized_model = "output/model_output_dir/#{File.basename(patient.step(:personalize_patient_BB).inputs[:model_prefix])}_personalized"
      options[:bnd_file] = patient.file(personalized_model + '.bnd')
      options[:cfg_file] = patient.file(personalized_model + '.cfg')
      jobs = [{:inputs => options.dup}]

      maboss.produce.file('output').ko_file.list.each do |ko|
        personalized_model = "output/model_output_dir/#{File.basename(patient.recursive_inputs[:model_prefix])}_personalized__#{ko}_ko"
        options[:bnd_file] = patient.file(personalized_model + '.bnd')
        options[:cfg_file] = patient.file(personalized_model + '.cfg')
        jobs << {:inputs => options.dup, :jobname => jobname + "_#{ko}_ko"}
      end

      jobs
    end
  task :PhysiBoSS => :array do
    dependencies.select{|d| d.task_name == :PhysiBoSS_BB}.collect{|d| d.load }.flatten
  end


  desc <<~EOF
  Conduct a meta-analysis over a set of patients
  EOF
  input :meta_file, :path, "Metadata file", nil, :nofile => true
  input :repetitions, :integer, "Repetitions for PhysiBoSS", 5
  dep :PhysiBoSS, :max_time => 100, :p_group => :placeholder, :repetition => :placeholder do |jobname,options|

    metadata_file = options[:meta_file]
    if TSV === metadata_file
      tsv = metadata_file
    elsif (Step === metadata_file && metadata_file.load.start_with?("#"))
      tsv = metadata_file.produce.path.tsv :type => :list
    else
      tsv = TSV.open(metadata_file, :type => :list)
    end

    jobs = tsv.collect do |id,values|
      group, file = values
      options[:repetitions].times.collect do |rep|
        {:inputs => options.merge(:p_group => group, :p_file => Path.setup(file), :repetition => rep), :jobname => id, :id => id}
      end
    end.flatten.compact

    jobs
  end
  task :meta_analysis => :array do 

    results_dir = file('results_dir')
    rec_dependencies
      .select{|d| d.task_name == :PhysiBoSS_BB }
      .each do |phy_bb|
        sample = phy_bb.inputs[:sample]
        model = File.basename(phy_bb.inputs[:bnd_file]).sub('.bnd','')
        rep = phy_bb.recursive_inputs[:repetition]
        filename = [model, 'physiboss_run', rep.to_i + 1] * "_"
        target = results_dir[sample].physiboss_results[filename]
        Open.cp phy_bb.file('output').results_dir, target
      end

    options = inputs.to_hash
    meta_file = file('metadata.tsv')
    Open.write(meta_file, options[:meta_file].to_s(header_hash: '', preamble: false))
    options[:meta_file] = meta_file
    options[:ko_file] = step(:MaBoSS_BB).join.file('output').ko_file
    options[:reps] = options[:repetitions]
    options[:model_prefix] = recursive_inputs[:model_prefix]
    options[:out_dir] = results_dir
    options[:verbose] = 1
    job = PerMedCoE.job(:meta_analysis_BB, clean_name, options)
    job.clean if job.error?
    job.produce
    job.file('output').glob("**/*")
  end

  dep :sample_info
  dep_task :pilot, Covid19, :meta_analysis, :meta_file => :sample_info
end

require 'rbbt/tasks/Covid19'
