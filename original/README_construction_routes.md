# Two strongly regular graphs with parameters (768, 312, 136, 120)

18 November 2026. Both constructions realize the complementary pair that Brouwer's public table marked open when checked on this date. The Cayley example gives a regular partial difference set in a nonabelian group of order 768. The second graph has exactly three vertex orbits. [Brouwer's table](https://aeb.win.tue.nl/graphs/srg/srgtab751-800.html)

| Witness | Four-cliques through a vertex in each of the three fibres | SHA-256 of graph6 |
|---|---|---|
| [Cayley graph](../graphs/srg768_312_136_120_cayley.g6) | 445,040 in every fibre | b22543c3d66965204f5b1b5100e1c439257e1ebaa891bdb01568687ee1c3b547 |
| [Three-orbit graph](../graphs/srg768_312_136_120_threeorbit.g6) | 444,230 / 443,250 / 443,760 | 267c53eebd256dd63a2a3deee7395b55a865ab68b844e1f00dc6a566c25f19b7 |
| [Complement of the three-orbit graph](../graphs/srg768_455_262_280_threeorbit_complement.g6) | Parameters $(768,455,262,280)$ | bb147cb5850ffb5138da2c595183b80e9d948136bcb4bdb4ad14a679254b33ed |

Translations act transitively within each fibre. The three different four-clique counts therefore distinguish the fibres of the second graph under every automorphism and prove that the two witnesses are nonisomorphic. The reconstruction counts triangles in each representative's neighbourhood. It also checks degrees, symmetry and common neighbours for every unordered pair. An independent graph6 decoder (`../code/verify_srg.py` in this repository) verified all three exports. [Reconstruction](construct768.py)

## The signed Walsh transform

Use $F=\mathbb F_2[t]/(t^4+t+1)$ and $V=F^2$. Integers from 0 through 15 encode field elements by the coefficients of $1,t,t^2,t^3$, in increasing bit order. Identify $V$ with its character group by $\chi_{(a,b)}(x,y)=(-1)^{\operatorname{Tr}_{F/\mathbb F_2}(ax+by)}$. Write $Hf(x)=\sum_{\beta\in V}f(\beta)\chi_\beta(x)$, so $H^2=256I$.

For $m$ in the ordered tuple $(0,1,2,10)$ take the 15 nonzero elements of the line $(a,ma)$. Partition these by multiplying $a$ by the order-five element encoded by 8. Within each line order the three orbits by the smallest integer representing their first coordinate. Number the resulting twelve orbits from 0 to 11. Define $z=-1$ on orbits $0,1,3,6,11$, and $z=1$ on the other seven. Set $z=0$ off these four punctured lines.

The signed sum defining $Hz$ is 10 at 96 positions and $-6$ at 160 positions. This finite identity can be checked in 256 sums of 60 signs with the definition above. Consequently $b=(Hz+6)/16$ is a Boolean function of weight 96 and $Hb=16z+96\delta_0$. This gives a cross-fibre adjacency block whose 60 nonconstant Fourier coefficients have absolute value 16. The search checked all $\binom{12}{5}=792$ sign choices on representatives of both semilinear orbits of four spread lines. Twelve choices worked on the representative used here. [Seed search](ports768.py) · [Full list of sign choices](ports768.jsonl)

## Assembly and spectrum

Index three copies of $V$ by $i\in\mathbb Z/3$. Label edge $\{i,i+1\}$ with a signed spectrum $z_i$. For the three-orbit graph apply the shear $(a,b)\mapsto(a,b+ha)$ to the seed for $h=0,4,13$, respectively, transporting its signs. The resulting four-line supports $S_i$ are pairwise disjoint. The same assembly works with the cyclic choice of supports described below.

For edge $i$, connect $x$ in fibre $i$ to $y$ in fibre $i+1$ when $b_i(x+y)=1$, where $b_i=(Hz_i+6)/16$. At fibre $i$ let $T_i=S_{i-1}\cup S_i$. This consists of eight punctured spread lines. Within that fibre use the connection set formed by the eight punctured annihilator lines of those eight lines. It has size 120. Its Fourier transform is $+8$ on $T_i$ and $-8$ at all other nonzero frequencies. Indeed the transform of a single punctured annihilator line is 15 on its dual line and $-1$ elsewhere.

Every vertex has $120+96+96=312$ neighbours. On the three-dimensional space of fibre constants the quotient has diagonal 120 and off-diagonal 96, hence eigenvalues $312,24,24$. Each nonzero frequency belongs to at most one $S_i$. If it belongs to none, its three-dimensional block is $-8I$. If it belongs to $S_i$, one diagonal entry is $-8$ and the remaining two-by-two block has diagonal 8 and off-diagonal $16z_i(\beta)$. Its eigenvalues are 24 and $-8$. The full spectrum is therefore $312^1,24^{182},(-8)^{585}$, giving

$$A^2=192I+16A+120J.$$

The graph has 136 common neighbours on each of its 119,808 adjacent pairs and 120 on each of its 174,720 nonadjacent pairs.

## The Cayley construction

On frequencies use $M(a,b)=(a,5a+6b)$ in the field encoding above. This map has order three. Transport the seed along $I,M,M^2$. The three slope sets are $\{0,1,2,10\}$, $\{3,5,9,12\}$ and $\{6,8,11,15\}$, so the spectral proof applies. The inverse transpose $R(x,y)=(x+5(7y),7y)$ acts on vertices. Applying $R$ and cycling the fibres preserves adjacency. With translations of $V$, this gives the regular action of $V\rtimes_R C_3$.

In the implementation a group element is $(i,x)$ with product $(i,x)(j,y)=(i+j,x+R^i y)$. The neighbours $D$ of the identity form an inverse-closed set of size 312. The program regenerates every adjacency row from its left translate and counts all ordered quotients $de^{-1}$ for $d,e\in D$. The counts are 312 at the identity, 136 on $D$ and 120 on all remaining elements. The linear map $R$ has distinct eigenvalues 1 and 7 over $F$. Its fixed additive factor has size 16 and its other factor is acted on by an order-three scalar. GAP identifies the resulting direct product as SmallGroup(768,1085323). [Cayley verifier](construct768.py)

## Reproduction

Run from this directory. The graph construction and pair checks use Python 3.10 or later and the standard library. The optional search for the seed uses NumPy.

~~~sh
python -B construct768.py --output /tmp/rebuild-srg768
python -B ../code/verify_srg.py 768,312,136,120 /tmp/rebuild-srg768/srg768_312_136_120*.g6
python -B ports768.py
~~~

