module ExTRI2

  task tf_types: :tsv do 
    tsv = Rbbt.data['tmp.tf.type.tsv'].tsv key_field: 'human_symbol', :header_hash => '', fields: ['TF type'], type: :double, sep2: ';'
    tsv = tsv.to_single
    Log.tsv tsv
    tsv2 = Rbbt.data['tmp.tf.type.tsv'].tsv key_field: 'Symbol', :header_hash => '', fields: ['TF type'], type: :single
    tsv2.each do |tf,type|
      next unless tsv.include? tf
      tsv[tf] = type if type
    end
    tsv
  end

  dep :CollecTRI2
  dep :tf_types
  task :CollecTRI2_tf_types => :tsv do
    tsv = step(:CollecTRI2).load
    types = step(:tf_types).load

    tsv.add_field 'TF Type' do |pair,values|
      tf = pair.split(':').first
      types[tf]
    end


    tsv
  end

  task :test => :yaml do
  end
end
