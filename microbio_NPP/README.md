# 1. Context

On souhaite estimer le niveau de contamination, par une bactérie
pathogène donnée, d’un produit alimentaire donné, au sein de la
production française. Pour cela on tire au hasard 50 produits issus de
la production française. Chaque produit fait l’objet d’une analyse
microbiologique suivant la méthode dite du Nombre le Plus Probable
(NPP). Cette méthode consiste à mettre en culture des échantillons du
produit dilués à trois dilutions décimales successives. En pratique,
pour chaque produit, sont mis en culture r = 3 tubes contenant m1 = 1g
du produit, r = 3 tubes contenant m2 = 0.1g du produit et r = 3 tubes
contenant m3 = 0.01g du produit. Le résultat brut de l’analyse
microbiologique effectuée sur chaque produit est constitué par les
nombres de tubes positifs à chacune des trois dilutions : nplus1, nplus2
et nplus3. Pour analyser les résultats pour chaque produit, c’est-à-dire
estimer la concentration en bactéries dans le produit lorsque c’est
possible, on fera l’hypothèse que les bactéries sont réparties de façon
homogène au départ dans le produit. On peut donc caractériser le nombre
de bactéries n pour une masse m prélevée dans un produit caractérisé par
une concentration en bactéries conc, par la loi de Poisson:
$P (n=x) = e^{-m \\times conc} \\times \\frac{(m \\times conc)^x}{x!}$

# 2. Model formalization

![dag](./assets/Disegno%20senza%20titol_microbio_npp.png)
