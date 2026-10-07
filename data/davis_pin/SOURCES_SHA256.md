# DAVIS raw inputs (DeepDTA mirror), fetched 2026-10-07
Base: https://raw.githubusercontent.com/hkmztrk/DeepDTA/master/data/davis/
Y is a pickled numpy array, shape (68, 442) = 68 ligands x 442 proteins (matches DAVIS size). Not yet cross-checked against our earlier DAVIS numbers.
9855c4f234ec6dd18295b0e0cd6ee54214cafa908c6eb1d447b47d454fc0294a  Y
141094d99510a52091d2ea1b05be7c1f66b2c5fcaeeeb261eaf4cbf2128bde6d  ligands_can.txt
187cedeee10cd3b58af915a2861a5ed0c4ff42a2d728f4dea5e0b47e3d75b291  proteins.txt
bb7a63b2178c40a4e2417bc701b95511407787fa5eabb6a43a695fe6bd722f56  folds/train_fold_setting1.txt
c97589fa2fae318b01870854af0ab6fdd155dd38524cc6ac1d0319418ac6d568  folds/test_fold_setting1.txt
Status: sources pinned only. No collision floor computed. R3 stays closed negative.
Family labels (H4): http://www.kinhub.org/kinases.html fetched 2026-10-08, saved as data/davis_pin/kinhub_kinases.html
c81227b4ac53c732b508f3c4f1d3d55fe491ac3b873516f7f4eeea42dc09b1c6  kinhub_kinases.html  (columns: xName, Manning name, HGNC, kinase name, Group, Family, SubFamily, UniProt)
