require 'rbbt/vector/model/huggingface'
require 'rbbt/vector/model/pytorch_lightning'
module ExTRI2

  input :base_model, :string, "Base model to start from", 'microsoft/BiomedNLP-BiomedBERT-large-uncased-abstract'
  input :folds, :integer, "Cross-validation folds", 1
  task :train_TRI_tag => :tsv do |base_model,folds|
    data = Rbbt.classifiers_training.data["TRI_data.tsv"].find

    tokenizer_options = 
      { 
        truncation: true,
        padding: true,
        model_max_length: 512,
      }

    model_options = { }

    training_options = 
      {
        learning_rate: 5e-6,
        num_train_epochs: 3,
        save_strategy: 'no',
        report_to: "none",
      }

    model = HuggingfaceModel.new "SequenceClassification", base_model, files_dir, 
      training_options: training_options, tokenizer_options: tokenizer_options, model_options: model_options

    model.extract_features do |sentence,tf,tg|
      sentence.gsub("[TF]", "<TF>#{tf}</TF>").gsub("[TG]", "<TG>#{tg}</TG>")
    end

    traverse data, type: :array do |line|
      id, sentence, tf, tg, valid, mor, validated = line.split("\t")[2..-1]
      next unless id.include?(":")
      model.add [sentence, tf, tg], valid.to_i
    end

    log :training, "Starting cross-validation with #{folds} folds"
    model.cross_validation(folds)
  end

  input :base_model, :string, "Base model to start from", 'microsoft/BiomedNLP-BiomedBERT-large-uncased-abstract'
  input :folds, :integer, "Cross-validation folds", 1
  task :train_TRI_masked => :tsv do |base_model,folds|
    data = Rbbt.classifiers_training.data["TRI_data.tsv"].find

    tokenizer_options = 
      { 
        truncation: true,
        padding: true,
        model_max_length: 512,
      }

    model_options = { }

    training_options = 
      {
        learning_rate: 5e-6,
        num_train_epochs: 3,
        save_strategy: 'no',
        report_to: "none",
      }

    model = HuggingfaceModel.new "SequenceClassification", base_model, files_dir, 
      training_options: training_options, tokenizer_options: tokenizer_options, model_options: model_options

    #model.extract_features do |sentence,tf,tg|
    #  sentence.gsub("[TF]", "<TF>#{tf}</TF>").gsub("[TG]", "<TG>#{tg}</TG>")
    #end

    traverse data, type: :array do |line|
      id, sentence, tf, tg, valid, mor, validated = line.split("\t")[2..-1]
      next unless id.include?(":")
      model.add sentence, valid.to_i
    end

    log :training, "Starting cross-validation with #{folds} folds"
    model.cross_validation(folds)
  end

  input :base_model, :string, "Base model to start from", 'microsoft/BiomedNLP-BiomedBERT-large-uncased-abstract'
  input :folds, :integer, "Cross-validation folds", 1
  input :learning_rate, :float, "Learning rate", 5e-6
  input :epochs, :integer, "Number of epochs", 3
  task :train_TRI => :tsv do |base_model,folds,learning_rate,epochs|
    data = Rbbt.classifiers_training.data["TRI_data.tsv"].find

    tokenizer_options = 
      { 
        truncation: true,
        padding: true,
        model_max_length: 512,
      }

    model_options = { }

    training_options = 
      {
        learning_rate: learning_rate,
        num_train_epochs: epochs,
        save_strategy: 'no',
        report_to: "none",
      }

    model = HuggingfaceModel.new "SequenceClassification", base_model, files_dir, 
      training_options: training_options, tokenizer_options: tokenizer_options, model_options: model_options

    m, tok = model.init

    new_tokens = ["[TF]", "[TG]"]
    tok.add_tokens(new_tokens)
    m.resize_token_embeddings(PyCall.len(tok))

    traverse data, type: :array do |line|
      id, sentence, tf, tg, valid, mor, validated = line.split("\t")[2..-1]
      next unless id.include?(":")
      model.add sentence, valid.to_i
    end

    log :training, "Starting cross-validation with #{folds} folds"
    model.cross_validation(folds)
  end

  input :base_model, :string, "Base model to start from", 'microsoft/BiomedNLP-BiomedBERT-large-uncased-abstract'
  input :folds, :integer, "Cross-validation folds", 1
  input :learning_rate, :float, "Learning rate", 5e-6
  input :epochs, :integer, "Number of epochs", 3
  task :train_custom => :tsv do |base_model,folds,learning_rate,epochs|
    data = Rbbt.classifiers_training.data["TRI_data.tsv"].find

    tokenizer_options = 
      { 
        model_max_length: 150,
      }

    model_options = { }

    training_options = 
      {
        batch_size: 64,
        max_epochs: 1,
      }

    tok = RbbtPython.import_method("rbbt_dm.huggingface", :load_tokenizer).call(base_model, **tokenizer_options)
    new_tokens = ["[TF]", "[TG]"]
    dict = {'additional_special_tokens' => new_tokens}
    tok.add_special_tokens(special_tokens_dict: dict)

    vocab_size = PyCall.len(tok)

    model = PytorchLightningModel.new files_dir, "TRICustomModel", "tri_token_classification",
      batch_size: 1,
      training_options: training_options,
      pretrained_model_name: base_model, 
      tf_token_id: vocab_size - 2,
      tg_token_id: vocab_size - 1,
      vocab_size: vocab_size

    model.extract_features do |sentence,sentences|
      emb = if sentences
              tok.call(sentences, return_tensors: 'pt', padding: true, truncation: true)
            else
              tok.call(sentence, return_tensors: 'pt', padding: true, truncation: true)
            end
      features = RbbtPython.numpy2ruby(emb["input_ids"]).zip(RbbtPython.numpy2ruby(emb["attention_mask"])).collect{|i,a| 
        i + a
      }
    end

    model.post_process do |predictions|
      RbbtPython.collect(predictions) do |logits|
        logits = RbbtPython.numpy2ruby logits
        best_class = logits.index logits.max
        best_class = model_options[:class_labels][best_class] if model_options[:class_labels]
        best_class
      end
    end

    sentences = []
    labels = []
    i = 100
    traverse data, type: :array do |line|
      id, sentence, tf, tg, valid, mor, validated = line.split("\t")[2..-1]
      next if sentence.nil? || sentence.strip.empty?
      parts = sentence.split(/[^a-z\[\]\(\)0-9]/i)
      next if parts.length > 100
      next unless parts.include?("[TF]")
      next unless parts.include?("[TG]")
      next unless id.include?(":")
      label = valid.to_i
      sentences << sentence
      labels << label
      i -= 1
      break if i == 0
    end
    model.add_list sentences, labels

    log :training, "Starting cross-validation with #{folds} folds"
    model.cross_validation(folds)
  end
end
