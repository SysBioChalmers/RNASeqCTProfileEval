%export metabolic genes from HumanGEM
HumanGEMPath = HumanGemInstall.getHumanGEMPath();
load(strcat(HumanGEMPath, '/ModelFiles/mat/HumanGEM.mat'))
[grRules,genes,rxnGeneMat] = translateGrRules(ihuman.grRules, 'Name');
writecell(genes,'data/metabolic_genes.txt');