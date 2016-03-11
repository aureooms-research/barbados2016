import sys
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter

G = nx.Graph()

n = int( sys.argv[1] )
r = n

# n**2 vertices per round
for i in range( r ) :
    for x in range( n ) :
        for y in range( n ) :
            G.add_node( ( i , ( x , y ) ) )

# add increasing constraint
for i in range( r ) :
    for x1 in range( n ) :
        for y1 in range( n ) :
            for x2 in range( x1 , n ) :
                for y2 in range( y1 , -1 , -1 ) :
                    if x1 == x2 and y1 == y2 : continue
                    G.add_edge( ( i , ( x1 , y1 ) ) , ( i , ( x2 , y2 ) ) )
            for x2 in range( x1 , -1 , -1 ) :
                for y2 in range( y1 , n ) :
                    if x1 == x2 and y1 == y2 : continue
                    G.add_edge( ( i , ( x1 , y1 ) ) , ( i , ( x2 , y2 ) ) )

# add killing constraint
for i in range( r ) :
    for x1 in range( n ) :
        for y1 in range( n ) :
            for j in range( i + 1 , r ) :
                G.add_edge( ( i , ( x1 , y1 ) ) , ( j , ( x1 , y1 ) ) )
                for x2 in range( x1 + 1 , n ) :
                    G.add_edge( ( i , ( x1 , y1 ) ) , ( j , ( x2 , y1 ) ) )
                for y2 in range( y1 + 1 , n ) :
                    G.add_edge( ( i , ( x1 , y1 ) ) , ( j , ( x1 , y2 ) ) )

degree_sequence=sorted(nx.degree(G).values()) # degree sequence
#print "Degree sequence", degree_sequence
dmax=max(degree_sequence)
dmin=min(degree_sequence)
davg=sum(degree_sequence)/len(degree_sequence)
dmed=degree_sequence[len(degree_sequence)//2]
p32=int(n**(3/2))-1
d32=degree_sequence[p32]
print( 'number of vertices: {}'.format( len(degree_sequence) ) )
print( 'minimum degree: {}'.format( dmin ) )
print( 'median degree: {}'.format( dmed ) )
print( 'maximum degree: {}'.format( dmax ) )
print( 'average degree: {}'.format( davg ) )
print( '"n^(3/2)"th degree: {}'.format( d32 ) )
print( '"n^(3/2)" first degrees sum: {}'.format( sum( degree_sequence[:p32+1] ) ) )

if False :

    degree_counter=Counter(nx.degree(G).values()) # degree sequence
    # degree_sequence=sorted(nx.degree(G).values(),reverse=True) # degree sequence
    #print "Degree sequence", degree_sequence
    dmax=max(degree_sequence)

    # plt.loglog(degree_sequence,'b-',marker='o')
    plt.plot(*zip(*sorted(degree_counter.items())),'bo')
    plt.title("Degree histogram")
    plt.ylabel("count")
    plt.xlabel("degree")
    plt.plot([d32],[degree_counter[d32]],'ro')
    plt.plot([dmin],[degree_counter[dmin]],'go')
    plt.plot([dmed],[degree_counter[dmed]],'go')
    plt.plot([dmax],[degree_counter[dmax]],'go')

    # draw graph in inset
    # plt.axes([0.45,0.45,0.45,0.45])
    # Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]
    # pos=nx.spring_layout(Gcc)
    # plt.axis('off')
    # nx.draw_networkx_nodes(Gcc,pos,node_size=20)
    # nx.draw_networkx_edges(Gcc,pos,alpha=0.4)

    plt.savefig("degree-histogram-{}.svg".format(n))
    plt.show()

# take complement
H = nx.complement( G )

# find cliques
cliques = nx.find_cliques( H )
score = max( map( len , cliques ) )
print( 'best solutions score : {}'.format( score ) )

cliques = nx.find_cliques( H )
for i , clique in enumerate( filter(lambda c : len( c ) == score , cliques ) , 1 ) :
    print( 'solution #{}'.format( i ) )
    for round , point in  sorted( clique ) :
        print( 'round {}, play point {}'.format( round , point ) )

    grid = zip( *( ( i // n , i % n ) for i in range( n**2 ) ) )
    plt.plot( *grid,'bo')
    points = zip( *( point for _, point in clique ) )
    plt.plot(*points, 'ro')
    for round , point in clique :
        plt.annotate(str(round), xy=point, textcoords='data')
    plt.axis([-10, n+9, -10, n+9])
    plt.savefig("sol-networkx/{}#{}.svg".format(n,hash(tuple(map(tuple,clique)))))
    plt.clf()
    # plt.show()
