import sys
import matplotlib.pyplot as plt

n = int( sys.argv[1] )
r = n

p = MixedIntegerLinearProgram(maximization=True, solver = "GLPK")
x = p.new_variable(binary=True)
p.set_objective(sum(x[i,j,k] for i in range(r) for j in range(n) for k in range(n)))

# add increasing constraint
for i in range( r ) :
    for x1 in range( n ) :
        for y1 in range( n ) :
            for x2 in range( x1 , n ) :
                for y2 in range( y1 , -1 , -1 ) :
                    if x1 == x2 and y1 == y2 : continue
                    p.add_constraint( x[i , x1 , y1] + x[i , x2 , y2] <= 1 )
            for x2 in range( x1 , -1 , -1 ) :
                for y2 in range( y1 , n ) :
                    if x1 == x2 and y1 == y2 : continue
                    p.add_constraint( x[i , x1 , y1] + x[i , x2 , y2] <= 1 )

# add killing constraint
for i in range( r ) :
    for x1 in range( n ) :
        for y1 in range( n ) :
            for j in range( i + 1 , r ) :
                p.add_constraint( x[i , x1 , y1 ] + x[ j , x1 , y1 ] <= 1 )
                for x2 in range( x1 + 1 , n ) :
                    p.add_constraint( x[i, x1, y1] + x[j, x2, y1 ] <= 1 )
                for y2 in range( y1 + 1 , n ) :
                    p.add_constraint( x[i, x1, y1] + x[j, x1, y2 ] <= 1 )

# solve MIP
p.show()
print 'Objective Value:', p.solve()
for i, v in p.get_values(x).iteritems():
    print 'x_%s = %s' % (i, int(round(v)))

# cliques = nx.find_cliques( H )
# score = max( map( len , cliques ) )
# print( 'best solutions score : {}'.format( score ) )

# cliques = nx.find_cliques( H )
# for i , clique in enumerate( filter(lambda c : len( c ) == score , cliques ) , 1 ) :
    # print( 'solution #{}'.format( i ) )
    # for round , point in  sorted( clique ) :
        # print( 'round {}, play point {}'.format( round , point ) )

    # grid = zip( *( ( i // n , i % n ) for i in range( n**2 ) ) )
    # plt.plot( *grid,'bo')
    # points = zip( *( point for _, point in clique ) )
    # plt.plot(*points, 'ro')
    # for round , point in clique :
        # plt.annotate(str(round), xy=point, textcoords='data')
    # plt.axis([-10, n+9, -10, n+9])
    # plt.savefig("sol/{}#{}.svg".format(n,hash(tuple(map(tuple,clique)))))
    # plt.clf()
    # plt.show()
