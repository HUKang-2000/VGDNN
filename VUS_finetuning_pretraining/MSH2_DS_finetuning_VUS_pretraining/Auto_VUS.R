library(tensorflow)
library(keras)

args <- commandArgs(trailingOnly = TRUE)

fold        <- as.numeric(args[1])
auto_noise  <- as.numeric(args[2])
auto_drop   <- as.numeric(args[3])
hidden_node <- as.numeric(args[4])
gpu         <- as.numeric(args[5])

Sys.setenv("CUDA_VISIBLE_DEVICES" = gpu)

load(file.path("Data",
               paste0(fold, "-th_Transfer_setting_MSH2_DiseaseSpecific_forNN_30sets_Nafill_standardized.RData")))

d<-paste("x_pre_train_",fold,"_MSH2",sep="")
x_train<-get(d)
x_train<-apply(x_train[,-1],2,as.numeric)

d<-paste("x_test_",fold,"_MSH2",sep="")
x_test<-get(d)
x_test<-apply(x_test[,-1],2,as.numeric)
#============================================================================================================
#============================================================================================================
#pretraining
set.seed(24)
input_layer<-layer_input(shape = c(length(colnames(x_train))))

Hidden1<-layer_dense(units = hidden_node)
Hidden2<-layer_dense(units = hidden_node)
Hidden3<-layer_dense(units = hidden_node)
Hidden4<-layer_dense(units = hidden_node)
Hidden5<-layer_dense(units = hidden_node)

Regress<-input_layer%>%
  layer_gaussian_noise(stddev = auto_noise)%>%
  Hidden1%>%
  layer_activation_leaky_relu(alpha=0.2)%>%
  layer_dropout(rate=auto_drop)%>%
  Hidden2%>%
  layer_activation_leaky_relu(alpha=0.2)%>%
  layer_dropout(rate=auto_drop)%>%
  Hidden3%>%
  layer_activation_leaky_relu(alpha=0.2)%>%
  layer_dropout(rate=auto_drop)%>%
  Hidden4%>%
  layer_activation_leaky_relu(alpha=0.2)%>%
  layer_dropout(rate=auto_drop)%>%
  Hidden5%>%
  layer_activation_leaky_relu(alpha=0.2)%>%
  layer_dropout(rate=auto_drop)

Regression_model<-keras_model(
  inputs =input_layer,
  outputs = layer_dense(units=length(x_train[1,]),activation='linear')(Regress)
)

Regression_model %>% compile(
  loss='mean_absolute_error',
  optimizer_adam(lr = 1e-4,decay = 1e-3)
)

fit_model<-fit(Regression_model,
               x_train,
               x_train,
               epochs=40,
               shuffle=TRUE,
               batch_size=200,
               validation_data = list(x_test,x_test)
)

save_model_weights_tf(
  Regression_model,
  filepath = file.path("VUS_finetuning_pretraining", "MSH2_DS_finetuning_VUS_pretraining", "pretraining_VUS",
                       paste0(fold, "-th_MSH2_pretraining_VUS_",auto_drop, "-", auto_noise,"_autoencoder.ckpt")
  )
)

q()