require 'rbbt/util/python'
module Covid19

  input :patient_models, :array, "Personalized patient models from single cell analysis"
  input :repetitions, :integer, "Repetitions for PhysiBoSS", 2
  task :pycompss_physiboss => :array do |patient_models, repetitions|

    pycompss_script =<<-EOF
from PhysiBoSS_BB import physiboss_model, physiboss
from pycompss.api.api import compss_wait_on_directory

def main():
  physiboss_results = []
    EOF

    max_time = 100
    patient_model_info = {}
    output_dir = file('results_dir')
    patient_models.each do |model_dir|
      Path.setup(model_dir)
      model_dir = model_dir.find
      sample = Misc.digest(model_dir)
      model_files = model_dir.glob("*.bnd")
      model_files.each do |model_file|
        prefix = File.basename(model_file).sub('.bnd','')
        repetitions.times do |rep|
          rep = rep + 1
          physiboss_results = output_dir[sample].physiboss_results[prefix + "_physiboss_run_#{rep}"]
          model_output_dir = output_dir[sample]
          Open.mkdir model_output_dir
          Open.mkdir model_output_dir
          Open.mkdir model_output_dir.tmp
          pycompss_script +=<<-EOF
# PHYSIBOSS
  results_dir = "#{physiboss_results}"
  physiboss_model(sample="#{sample}",
            repetition=#{rep},
            prefix="#{prefix}",
            model_dir="#{model_dir}",
            out_file="#{model_output_dir.output}",
            err_file="#{model_output_dir.error}",
            results_dir=results_dir,
            max_time=#{max_time},
            parallel=1,
            tmpdir="#{model_output_dir.tmp}")
  physiboss_results.append(results_dir)
          EOF
        end
      end
    end

    pycompss_script +=<<-EOF
  for physiboss_result in physiboss_results:
      compss_wait_on_directory(physiboss_result)

if __name__ == '__main__':
  main()
    EOF

    cpus = config(:cpus, :runcompss, :default => 4)
    ppp pycompss_script
    TmpFile.with_file(pycompss_script, :extension => 'py') do |script|
      CMD.cmd_log("runcompss #{script}")
    end
    output_dir.glob("*")
  end

  desc <<~EOF
  Conduct a meta-analysis over a set of patients
  EOF
  input :meta_file, :path, "Metadata file", nil, :nofile => true
  input :repetitions, :integer, "Repetitions for PhysiBoSS", 2
  dep :personalize_single_cell_patient, :p_group => :placeholder do |jobname,options|

    metadata_file = options[:meta_file]
    if Step === metadata_file || metadata_file.load.start_with?("#")
      tsv = metadata_file.produce.path.tsv :header_hash => '', :type => :list
    else
      tsv = TSV.open(metadata_file, :header_hash => '', :type => :list)
    end

    tsv.collect do |sample,values|
      group, file = values
      {:inputs => options.merge(:p_group => group, :p_file => file), :jobname => sample}
    end
  end
  dep :pycompss_physiboss do |jobname,options,dependencies|
    personalized = dependencies.flatten.select{|dep| dep.task_name.to_sym == :personalize_single_cell_patient }
    models = personalized.collect{|d| "" + d.file('output/model_output_dir') }
    {:inputs => options.merge(:patient_models => models), :jobname => jobname}
  end
  task :hybrid_meta_analysis => :array do 

    results_dir = file('results_dir')
    physi = dependencies.last
    personal = dependencies[0..-1]

    id2sample = {}
    personal.each do |dep|
      model_dir = dep.file('output').model_output_dir
      id = Misc.digest(model_dir)
      id2sample[id] = dep.clean_name
    end

    physi.file('results_dir').glob("*").each do |sample_dir|
      id = File.basename(sample_dir)
      sample = id2sample[id]
      sample_dir.physiboss_results.glob("*").each do |physi_dir|
        prefix = File.basename(physi_dir)
        target = results_dir[sample].physiboss_results[prefix]
        Open.cp physi_dir, target
      end
    end

    model_prefix = recursive_inputs[:model_prefix]
    model_prefix.find if Path === model_prefix

    options = inputs.to_hash
    options[:ko_file] = step(:MaBoSS_BB).join.file('output').ko_file
    options[:reps] = options[:repetitions]
    options[:model_prefix] = model_prefix
    options[:out_dir] = results_dir
    options[:verbose] = 0

    job = PerMedCoE.job(:meta_analysis_BB, clean_name, options)
    job.produce
    job.file('output').glob("**/*")
  end

  dep :sample_info
  dep_task :hybrid_pilot, Covid19, :hybrid_meta_analysis, :meta_file => :sample_info
end

