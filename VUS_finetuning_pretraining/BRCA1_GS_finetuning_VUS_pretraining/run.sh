GPU=0
auto_noise=1
auto_drop=0.7
hidden_node=9000
auto_epoch=20

for Fold in {1..30}
  do
    Rscript Auto_VUS.R $Fold $auto_noise $auto_drop $hidden_node $GPU
    for Ensemble in {1..50}
      do
        Rscript Auto_DS_Conti_Fine.R $Fold $Ensemble $auto_noise $auto_drop $auto_epoch $hidden_node $GPU
      done
  done