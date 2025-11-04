Workflow.require_workflow "SaezLab"

module ExTRI2
  input :source, :select, "Source of ExTRI2", :clean, select_options: %w(comnplete clean clean_no_MoR orig orig_sign extri1)
  dep :ExTRI2_clean, compute: :produce, only_authoritative_tfs: :placeholder, remove_auto_regulation: :placeholder do |jobname,options|
    case options[:source].to_sym
    when :complete
      {task: :ExTRI2_clean, inputs: options.merge(only_authoritative_tfs: false, remove_auto_regulation: false), jobname: jobname}
    when :clean
      {task: :ExTRI2_clean, inputs: options.merge(only_authoritative_tfs: false, remove_auto_regulation: true), jobname: jobname}
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
  task_alias :CollecTRI2_pre, ExTRI, :CollecTRI, not_overriden: true, confidence: 'none' do |jobname,options,dependencies|
    extri2 = dependencies.flatten.first
    options["ExTRI#ExTRI_final"] = extri2 if options[:source].to_s != "extri1"
    {jobname: jobname, inputs: options}
  end

  dep :CollecTRI2_pre
  task :CollecTRI2 => :tsv do
    tsv = step(:CollecTRI2_pre).load

    tsv.fields = tsv.fields.collect{|f| f.sub('ExTRI', 'ExTRI2')}

    tsv.process '[ExTRI2] present' do |v|
      v.first == 'ExTRI' ? ['ExTRI2'] : ''
    end

    tsv
  end
end
