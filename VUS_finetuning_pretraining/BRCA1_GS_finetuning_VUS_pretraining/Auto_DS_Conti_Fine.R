library(tensorflow)
library(keras)
library(PRROC)

args <- commandArgs(trailingOnly = TRUE)

fold        <- as.numeric(args[1])
model       <- as.numeric(args[2])
auto_noise  <- as.numeric(args[3])
auto_drop   <- as.numeric(args[4])
auto_epoch  <- as.numeric(args[5])
hidden_node <- as.numeric(args[6])
gpu         <- as.numeric(args[7])

Sys.setenv("CUDA_VISIBLE_DEVICES" = gpu)

load(file.path("Data",
               paste0(fold, "-th_Transfer_setting_BRCA1_DiseaseSpecific_forNN_30sets_Nafill_standardized.RData")))

d<-paste("x_train_",fold,"_BRCA1",sep="")
t<-paste("y_train_",fold,"_BRCA1",sep="")

x_train<-get(d)
x_train<-apply(x_train[,-1],2,as.numeric)
y_train<-get(t)

d<-paste("x_test_",fold,"_BRCA1",sep="")
t<-paste("y_test_",fold,"_BRCA1",sep="")

x_test<-get(d)
x_test<-apply(x_test[,-1],2,as.numeric)
y_test<-get(t)
#============================================================================================================
#============================================================================================================
#pretraining
set.seed(21)
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

load_model_weights_tf(
  Regression_model,
  filepath = file.path("VUS_finetuning_pretraining","BRCA1_GS_finetuning_VUS_pretraining", "pretraining_VUS",
                       paste0(fold, "-th_BRCA1_pretraining_VUS_",auto_drop, "-", auto_noise,"_autoencoder.ckpt")))

fit_model<-fit(Regression_model,
               x_train,
               x_train,
               epochs=auto_epoch,
               shuffle=TRUE,
               batch_size=100,
               validation_data = list(x_test,x_test)
)

d<-paste("x_train_",fold,"_BRCA1",sep="")
t<-paste("y_train_",fold,"_BRCA1",sep="")

x_train<-get(d)
BRCA1_train_idx<-which(x_train[,"train_GI"]=="BRCA1")
x_train<-apply(x_train[BRCA1_train_idx,-1],2,as.numeric)
y_train<-get(t)
y_train<-y_train[BRCA1_train_idx]

d<-paste("x_test_",fold,"_BRCA2",sep="")
t<-paste("y_test_",fold,"_BRCA2",sep="")

x_test<-get(d)
BRCA1_test_idx<-which(x_train[,"test_GI"]=="BRCA1")
x_test<-apply(x_test[BRCA1_test_idx,-1],2,as.numeric)
y_test<-get(t)
y_train<-y_train[BRCA1_test_idx]

#Fine_tuning
set.seed(21)
load(file.path("VUS_finetuning_pretraining","BRCA1_GS_finetuning_VUS_pretraining", "GS_Par_set_with_pretraining.RData"))

Dropout<-c(0.4,0.6,0.8)
Epoch<-c(10,20)
Batch<-c(10,20)

Parameter<-expand.grid(Dropout,Epoch,Batch)
colnames(Parameter)<-c("Dropout","Epochs","Batchs")

dropout<-Parameter[Par_set[fold],"Dropout"]
epochs<-Parameter[Par_set[fold],"Epochs"]
batchs<-Parameter[Par_set[fold],"Batchs"]

input_layer<-layer_input(shape = c(length(colnames(x_train))))

Classify<-input_layer%>%
  Hidden1%>%
  layer_activation_leaky_relu(alpha=0.2)%>%
  layer_dropout(rate=dropout)%>%
  Hidden2%>%
  layer_activation_leaky_relu(alpha=0.2)%>%
  layer_dropout(rate=dropout)%>%
  Hidden3%>%
  layer_activation_leaky_relu(alpha=0.2)%>%
  layer_dropout(rate=dropout)%>%
  Hidden4%>%
  layer_activation_leaky_relu(alpha=0.2)%>%
  layer_dropout(rate=dropout)%>%
  Hidden5%>%
  layer_activation_leaky_relu(alpha=0.2)%>%
  layer_dropout(rate=dropout)

Classification_model<-keras_model(inputs =input_layer,
                                  outputs = layer_dense(units=1,activation='sigmoid')(Classify)
)

Classification_model %>% compile(
  loss='binary_crossentropy',
  optimizer_adam(learning_rate = 1e-4,decay = 1e-2),
  metrics=metric_auc(curve = "PR",name = "AUPRC")
)

fit_model<-fit(Classification_model,
               x_train,
               y_train,
               epochs=epochs,
               shuffle=TRUE,
               batch_size=batchs,
               validation_data = list(x_test,y_test)
)

predictions_tr <- predict(Classification_model, as.matrix(x_train))
predictions_te <- predict(Classification_model,as.matrix(x_test))

fg_tr <- predictions_tr[y_train == 1]
bg_tr <- predictions_tr[y_train == 0]

fg_te <- predictions_te[y_test == 1]
bg_te <- predictions_te[y_test == 0]

Total_pr_tr <- pr.curve(scores.class0 = fg_tr, scores.class1 = bg_tr, curve = T)$auc.integral
Total_pr_te <- pr.curve(scores.class0 = fg_te, scores.class1 = bg_te, curve = T)$auc.integral

predicted_te<-predictions_te
predicted_tr<-predictions_tr

save(Total_pr_te, Total_pr_tr, predicted_tr, predicted_te,
     file = file.path("VUS_finetuning_pretraining","BRCA1_GS_finetuning_VUS_pretraining", "result", 
                      paste0("BRCA1_GS_", fold, "-", model, "-th_folds.RData")))

q()