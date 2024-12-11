require 'rbbt-util'
require 'rbbt/workflow'

Workflow.require_workflow "SaezLab"

module ExTRI2

  input :source, :select, "Source of ExTRI2", :clean, select_options: %w(clean clean_no_MoR orig orig_sign extri1)
  dep :ExTRI2_clean, compute: :produce do |jobname,options|
    case options[:source].to_sym
    when :complete
      {task: :ExTRI2_clean, inputs: options.merge(only_authoritative_tfs: false, remove_auto_regulation: false), jobname: jobname}
    when :clean
      {task: :ExTRI2_clean, inputs: options, jobname: jobname}
    when :clean_no_MoR
      {task: :ExTRI2_clean, inputs: options.merge(no_MoR: true), jobname: jobname}
    when :orig
      {task: :ExTRI2_clean_orig, inputs: options, jobname: jobname}
    when :orig_sign
      {task: :ExTRI2_clean_orig_sign, inputs: options, jobname: jobname}
    when :extri1
      nil
    end
  end
  task_alias :CollecTRI2, ExTRI, :CollecTRI, not_overriden: true do |jobname,options,dependencies|
      extri2 = dependencies.flatten.first
      options["ExTRI#ExTRI_final"] = extri2 if options[:source].to_s != "extri1"
      {jobname: jobname, inputs: options}
    end

  dep :CollecTRI2
  task_alias :regulome, SaezLab, :regulome, not_overriden: true do |jobname,options,dependencies|
    collectri2 = dependencies.flatten.first
    options["ExTRI#CollecTRI"] = collectri2 if options[:source].to_s != "extri1"
    {jobname: jobname, inputs: options}
  end

  dep :regulome
  task_alias :evaluate_knocktf, SaezLab, :evaluate_knocktf, "SaezLab#regulome" => :regulome, not_overriden: true

  dep :ExTRI2_clean
  dep SaezLab, :evaluate_knocktf, not_overriden: true do |jobname,options,dependencies|
    extri2 = dependencies.flatten.first
    jobs = []

    #[1,3,5,10].each do |support_evidence_max|
    #  [1,3,5,10].each do |sign_support_evidence_max|
    #    [0.1,0.25,0.5,0.75,0.9].each do |sign_support_balance|
    #      [0.1,0.25,0.5,0.75,0.9].each do |strict_negative|
    #        [5, 10, 20].each do |min_n|
    #          [true,false].each do |overriden|
    #            [:ulm, :mlm].each do |method|

    #[5,10,15].each do |support_evidence_max|
    #  [5,10,15].each do |sign_support_evidence_max|
    #    [0.1,0.5].each do |sign_support_balance|
    #      [0.1,0.5,0.75].each do |strict_negative|
    #        [5, 10, 20].each do |min_n|
    #          [false,true].each do |overriden|
    #            [:ulm, :mlm].each do |method|
     
    [0.1,0.5,0.9].each do |strict_negative|
      [0.1,0.5,0.9].each do |sign_support_balance|
        [10,50,100].each do |support_evidence_max|
          [0.1,0.5,0.9].each do |auto_regulation_weight|
            [3, 5, 10].each do |min_n|
              [true, false].each do |overriden|
                [:ulm, :mlm].each do |method|

                  inputs = options.merge(
                    strict_negative: strict_negative,
                    sign_support_balance: sign_support_balance,
                    support_evidence_max: support_evidence_max,
                    auto_regulation_weight: auto_regulation_weight,
                    min_n: min_n,
                    method: method
                  )

                  inputs["ExTRI#ExTRI_final"] = extri2 if overriden

                  jobs << { inputs: inputs, not_overriden: true }
                end
              end
            end
          end
        end
      end
    end
    jobs
  end
  task :knocktf_sweep => :tsv do
    dumper = TSV::Dumper.new(key_field: "Name", fields: "Missing,Miss-matches,Matches,Accuracy %,Average Rank,Version,strict_negative,sign_support_balance,support_evidence_max,auto_regulation_weight,min_n,method".split(","), type: :list, cast: :to_f)
    dumper.init
    TSV.traverse dependencies[1..-1], type: :array, bar: progress_bar, into: dumper do |dep|
      strict_negative = dep.recursive_inputs[:strict_negative]
      sign_support_balance = dep.recursive_inputs[:sign_support_balance]
      support_evidence_max = dep.recursive_inputs[:support_evidence_max]
      auto_regulation_weight = dep.recursive_inputs[:auto_regulation_weight]
      method = dep.recursive_inputs[:method]
      min_n = dep.recursive_inputs[:min_n]
      overriden = dep.step(:ExTRI_final).path.include?("ExTRI2")
      version = overriden ? "ExTRI2" : "ExTRI1"
      res = dep.load.values + [version, strict_negative,sign_support_balance,support_evidence_max,auto_regulation_weight,min_n,method]
      name = Misc.obj2digest(res)
      [name, res]
    end
  end


end
