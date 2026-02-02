# Projet - Implementation de McEliece en C (codes de Goppa)

**UE : Cryptographie post-quantique & theorie des codes**

Reference principale : `guide.pdf` (Sumi–Morozov–Takagi, 2011). Les algorithmes 1 a 7 de ce guide servent de fil conducteur.

## Objectif
Implementer un mini-cryptosysteme McEliece (codes de Goppa irreductibles) en C, de la generation de cle a la decryption.
Pour ce projet, **on ne met pas en oeuvre les matrices S et P** (pas de masquage ni permutation publique). On travaille donc directement avec G construit depuis H' systematisee.
Le projet est decoupe en sous-problemes mathematiques clairement identifiables. Les fichiers `.h` fournis listent **uniquement** les fonctions necessaires, avec les operations mathematiques attendues.

## Ce que vous devez produire
- Une implementation C propre des fonctions declarees dans les `.h`.
- Des tests simples pour valider chaque etape (voir section Tests).
- Un court rapport (2-3 pages) expliquant vos choix d'implementation et les difficultes.

## Organisation des fichiers
- `param.h` : definition des parametres (m, t, n, k) et polynome irreductible de F2^m.
- `gf.h` : arithmetique dans $F_{2^m}$ (addition, multiplication, inverse, racine carree).
- `poly.h` : polynomes sur $F_{2^m} (add/mul/mod/xgcd, evaluation, racine, etc.).
- `matrix.h` : matrices binaires $F_2$ (bitset), operations de lignes/colonnes.
- `key_gen.h` : generation des cles (Algorithme 1 du guide).
- `encrypt.h` : chiffrement (Algorithme 2).
- `decrypt.h` : decryption + decodage Patterson (Algorithmes 3 a 7).
- `util.h` : helpers (vecteurs binaires, poids de Hamming).
- `rng.h` : interface de tirage aleatoire.

## Description mathematique 

### Parametres
- m : degre du corps $F_2^m$ (n = 2^m)
- t : capacite de correction d'erreurs
- $k = n - m*t$
- $g(X)$ : polynome de Goppa irreductible de degre t sur $F_2^m$
- $L = (\alpha_0, ..., \alpha_{n-1})$ : support (elements distincts de $F_2^m$)

### Algorithme 1 (KeyGen) - a realiser dans `key_gen.h`
1. Construire la matrice de parite H = Y * Z, avec
   - $Y_{i,j} = \alpha_j^i$ pour $i = 0, \cdots, t-1$
   - $Z = diag(g(\alpha_j)^{-1})$
2. Developper H (coefficients dans $F_2^m$) en H' binaire de taille (m*t) x n.
3. Mettre H' en forme systematique $H^'_r = [A | I_{m*t}]$.
4. Construire $G = [I_k | A^T]$.
5. On **ne** tire **pas** S ni P pour ce projet.
6. $G_{pub} = G$.

### Algorithme 2 (Encrypt) - `encrypt.h`
- $y = m  G_{pub}$
- $c = y XOR e$, avec $w_H(e) = t$

### Algorithmes 3-7 (Decrypt + Decode) - `decrypt.h`
1. Decoder via Patterson (Algorithmes 4-7) pour trouver l'erreur e et y = c XOR e.
2. Recuperer m directement depuis G (pas de S, pas de P).

#### Patterson
Pour ce projet, vous devez **decomposer** l'algorithme en sous-fonctions claires (voir `decrypt.h`).  
On suit la structure des Algorithmes 4 a 7 du guide :

**Algorithme 4 : Decodage d'un code de Goppa binaire**
1. **Syndrome** :  
   $$S_c(X) = \\sum_{i=0}^{n-1} c_i (X - L_i)^{-1} \\bmod g(X).$$
2. **Test codeword** : si $$S_c(X)=0$$ alors $$y=c$$.
3. **Inverse** :  
   $$T(X) = S_c(X)^{-1} \\bmod g(X).$$
4. **Racine carree** :  
   $$\\theta(X) = \\sqrt{T(X) + X} \\bmod g(X)$$ (Algorithme 5).
5. **Equation cle** :  
   resoudre $$\\sigma(X)$$ via l'algorithme d'Euclide etendu (Algorithme 6).
6. **Recherche de racines** :  
   trouver les racines de $$\\sigma(X)$$ dans $$\\mathbb{F}_{2^m}$$ (Algorithme 7).
7. **Construction de l'erreur** :  
   $$e_i = 1$$ si $$L_i$$ est racine, sinon $$e_i = 0$$.
8. **Sortie** : $$y = c \\oplus e$$, sinon echec si aucune racine.

**Algorithme 5 : Racine carree modulo g(X)**  
Implementer le calcul de $$\\sqrt{Q(X)} \\bmod g(X)$$ avec pre-calculs (tables ou polynomes $$R_i$$).

**Algorithme 6 : Equation cle**  
Implementer l'algorithme d'Euclide etendu pour obtenir  
$$\\sigma(X) = \\gamma(X)^2 + X \\phi(X)^2.$$

**Algorithme 7 : Recherche de racines**  
Implementer un algorithme de type Cantor–Zassenhaus (ou variante) pour factoriser \( \\sigma(X) \) et extraire les racines.

## Tests minimaux
1. **GF** : verifie que $a * a^{-1} = 1$ pour $a \ne 0$.
2. **Poly** : verifie $p = (q * d + r)$ avec $deg(r) < deg(d)$.
3. **Systematisation** : verifie que le bloc droit de $H'_r$ est $I_{m*t}$.
4. **Orthogonalite** : verifie $H'_r * G^T = 0$.
5. **Chiffrement/Dechiffrement** : pour plusieurs messages, decrypt(encrypt(m)) = m.
6. **Patterson** : verifier que le decodeur corrige des erreurs de poids <= t sur des mots de code valides.

## Parametres conseilles pour les tests
Pour garder des temps raisonnables en TP, utilisez un petit jeu de parametres (ex. m = 6, t = 5), puis un jeu plus grand pour la validation finale.

## Bareme
- 25% : GF(2^m) + polynomes (operations de base correctes).
- 30% : KeyGen (H, H', systematisation, G).
- 25% : Chiffrement + Dechiffrement complet.
- 10% : Tests robustes et automatises.
- 10% : Rapport clair.

---

## Mapping fonctions -> operations mathematiques (resume)
Les commentaires de chaque fonction dans les `.h` indiquent precisement :
- les objets mathematiques manipules ($F2^m$, $F_2$, matrices, polynomes),
- la formule ou l'operation attendue,
- la reference a l'algorithme du guide (1 a 7).
