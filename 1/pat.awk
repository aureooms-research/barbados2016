#
# Exemple:
#
# ./gengraph pat 4 4 4 -check kindepsat 8 | ./glucose -model | awk -f pat.awk -v p=4 q=4 r=4
#
# SATISFIABLE
# p=4 q=4 r=4
#  .  2  .  1 
#  2  .  1  . 
#  .  .  4  3 
#  4  3  .  . 
# score=8
#

/^s/{
    model=substr($0,3);
    print model;
}

/^v/{ model=substr($0,3); }

END{
    if(match(model,/UNSAT/)) exit;

    print"p="p,"q="q,"r="r;
    split(model,V," "); # découpe le modèle en valeur
    N=p*q*r;

    for(i in V){
	v=V[i]; # v=variable=doit être dans [1,N]
	if((v<1)||(v>N)) continue;

	u=v-1;              # u=0..N-1=numéro de sommet
	z=int(u/(p*q));     # z=grille
	x=u%p;              # x=colonne
	y=int((u%(p*q))/p); # y=ligne

	G[x,y]=z+1;      # remplit la grille
	if(z>max) z=max; # calcule la valeur max
	score++;         # met à jour le score
    }
    
    # max=nombre de caractères pour écrire la valeur maximum de G
    max=length(""(max+1));

    for(y=q-1;y>=0;y--){
	for(x=0;x<p;x++)
	    printf(" %s ",((x,y) in G)? center(G[x,y],max) : ".");
	print"";
    }
    print"score="score;
    exit;
}

function center(x,t,
	     s,i,n){
# affiche x au centre d'une chaîne d'espaces de longueur t
    s=""x;
    n=t-length(s); # calcule le nombre d'espaces à mettre
    for(i=0;i<n;i++)
	if(i<n/2) s=" "s; # espace avant x
	else s=s" ";      # espace après x
    return s;
}
