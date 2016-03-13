import sys
import hashlib

# capital letter for a variable name means constant

M = int( sys.argv[1] ) # number of rows
N = int( sys.argv[2] ) # number of columns
R = int( sys.argv[3] ) # number of rounds

N = R * M * N # number of vertices of the dependence graph

lb = min( M , N ) # lower bound
ub = int( M * N ) # upper bound

best = None
nsol = 0

# k = round id
# i = row id
# j = column id

x = lambda k , i , j : 1 + k*(M*N) + i * N + j
Y00 = x(R-1,M-1,N-1) + 1
y = lambda i , b : Y00 + b*( N + 1 ) + i
# ijk = lambda v : ((v-1)//n**2,((v-1)%n**2)//n,(v-1)%n)

def out ( n , m , r , sol , mode = 'r' ) :

    h = hashlib.sha1(str(sol).encode('utf-8')).hexdigest()
    filepath = 'sol/{:0}x{}x{}/{}/{}'.format( m , n , r , val , h )

    try :
        os.makedirs( os.path.dirname( filepath ) )
    except :
        pass

    return open( filepath , mode )

def output_solution ( val , sol ) :

    global nsol
    nsol += 1


    print( 'solution #{} ({})'.format( nsol , path ) )
    print( 'solution score : {}'.format( val ) )

    moves = []
    for i in range(r):
        for j in range(n):
            for k in range(n):
                if sol[x(i,j,k)] :
                    moves.append( ( i + 1 , ( j , k ) ) )

    print( 'assert', val , '==', len(moves))
    # assert(val == len(moves))

    for round , point in moves :
        print( 'round {}, play point {}'.format( round , point ) )

    # plot
    grid = zip( *( ( i // n , i % n ) for i in range( n**2 ) ) )
    plt.plot(*(grid+['bo']))
    points = zip( *( point for _, point in moves ) )
    plt.plot(*(points+['ro']))
    for round , point in moves :
        plt.annotate(str(round), xy=point, textcoords='data')
    plt.axis([-10, n+9, -10, n+9])

    # save to file
    plt.savefig("{}.svg".format(path))

    # clear plot
    plt.clf()

    with open('{}.pat'.format(path),'w') as fd :
        fd.write('{} {}\n'.format(m , n , len(moves)))
        for round , point in moves :
            i, j = point
            fd.write('{} {} {}\n'.format(round,i,j))

def solve_puzzle ( pivot ) :

    p = SAT( solver = 'cryptominisat' )

    # add increasing constraint
    for i in range( r ) :
        for x1 in range( n ) :
            for y1 in range( n ) :
                for x2 in range( x1 , n ) :
                    for y2 in range( y1 , -1 , -1 ) :
                        if x1 == x2 and y1 == y2 : continue
                        p.add_clause( ( -x(i , x1 , y1) , -x(i , x2 , y2) ) )
                for x2 in range( x1 , -1 , -1 ) :
                    for y2 in range( y1 , n ) :
                        if x1 == x2 and y1 == y2 : continue
                        p.add_clause( ( -x(i , x1 , y1) , -x(i , x2 , y2) ) )

    # add killing constraint
    for i in range( r ) :
        for x1 in range( n ) :
            for y1 in range( n ) :
                for j in range( i + 1 , r ) :
                    p.add_clause( ( -x(i , x1 , y1) , -x( j , x1 , y1 ) ) )
                    for x2 in range( x1 + 1 , n ) :
                        p.add_clause( ( -x(i, x1, y1) , -x(j, x2, y1 ) ) )
                    for y2 in range( y1 + 1 , n ) :
                        p.add_clause( ( -x(i, x1, y1) , -x(j, x1, y2 ) ) )


    for i in range(1,N):
        for b in range(1,pivot+1):
            p.add_clause( ( -y( i , b ) , y( i - 1 , b ) , y( i-1, b-1) ) )
            p.add_clause( ( -y( i , b ) , y( i - 1 , b ) , i ) )

    for i in range(N):
        p.add_clause( ( y( i , 0 ) , ) )

    for b in range( 2 , pivot+1 ) :
        p.add_clause( ( -y( 0 , b ) , ) )

    p.add_clause( ( y( N , pivot ) , ) )

    return p()

while LB <= UB :

    pivot = ( LB + UB ) // 2

    print( 'LB {} UB {} pivot {} best {}'.format(LB , UB , pivot , best ))

    sol = solve_puzzle(pivot)

    if sol :
        output_solution(pivot,sol)
        best = pivot
        LB = pivot + 1
    else :
        UB = pivot - 1

