Each rda file containes three components:
gbm_models: the trained gradient boosting classifier
cutoff: less stringent threshold to be applied on Signature_3_mva column of SigMA outout. The pass_mva column is defined as the boolean: Signature_3_mva > cutoff.
cutoff_strict: stringent threshold to be applied on Signature_3_mva column of SigMA output. The pass_mva_strict is defined as the boolean: Signature_3_mva > cutoff_strict.
