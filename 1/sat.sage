import argparse
import hashlib

# capital letter for a variable name means constant


parser = argparse.ArgumentParser(description="Find Pat's game solutions for an M x N grid with R rounds.")
parser.add_argument('rows', metavar='M', type=int, help='number of rows')
parser.add_argument('columns', metavar='N', type=int, help='number of columns')
parser.add_argument('rounds', metavar='R', type=int, help='number of rounds')
parser.add_argument('-e', '--enumerate', action='store_true', help='enumerate all optimal solutions')
args = parser.parse_args()

M = args.rows # number of rows
N = args.columns # number of columns
R = args.rounds # number of rounds

best = None
best_sol = None
nsol = 0

# k = round id
# i = row id
# j = column id

x = lambda k , i , j : 1 + k*(M*N) + i * N + j
# ijk = lambda v : ((v-1)//n**2,((v-1)%n**2)//n,(v-1)%n)

def validate ( sol ) :
    pass

def out ( val , t , h , mode = 'r' ) :

    filepath = 'sol-web/{}x{}x{}/{}/{}/{}'.format( M , N , R , t , val , h )

    try :
        os.makedirs( os.path.dirname( filepath ) )
    except :
        pass

    print( 'opening "{}" with mode "{}"'.format( filepath , mode ) )

    return open( filepath , mode )

def getmoves ( sol ) :

    for k in range(R):
        for i in range(M):
            for j in range(N):
                if sol[x(k,i,j)] :
                    yield ( k + 1 , ( i + 1 , j + 1 ) )

def output_unsat ( val ) :

    sb = []
    sb.append('{} {} {}\n'.format(M , N , R))
    sb.append('{}\n'.format(val))
    sb.append('0\n') # UNSAT

    raw = ''.join(sb)
    h = hashlib.sha1(raw.encode('utf-8')).hexdigest()

    with out(val, 'unsat', h, 'w') as fd :
        fd.write(raw)

def output_solution ( val , sol ) :

    global nsol
    nsol += 1

    print( 'solution #{}'.format( nsol ) )
    print( '> score : {}'.format( val ) )

    moves = list(getmoves( sol ))

    if val != len(moves) :
        print( '> val {} != len(moves) {}'.format(val, len(moves)))
        print( '> still writing output just in case')

    for round , point in moves :
        print( '! round {}, play point {}'.format( round , point ) )


    sb = []
    sb.append('{} {} {}\n'.format(M , N , R))
    sb.append('{}\n'.format(val))
    sb.append('1\n') # SAT
    for round , point in moves :
        i, j = point
        sb.append('{} {} {}\n'.format(round,i,j))

    raw = ''.join(sb)
    h = hashlib.sha1(raw.encode('utf-8')).hexdigest()

    with out(val, 'sat', h, 'w') as fd :
        fd.write(raw)

def max_vars_sat ( p , nvars , firsty , variables , pivot ) :

    y = lambda k , t : firsty + t*( nvars + 1 ) + k

    # add 3-SAT clauses to maximize sum of x variables
    # y_{k,t} means sum_{i=1}^{k} x_i >= t

    # y_{0,0} is always true
    p.add_clause( ( y( 0 , 0 ) , ) )

    # y_{0,t} with t > 0 is always false
    for t in range(1, pivot+1 ):
        p.add_clause( ( -y( 0 , t ) , ) )

    for k , var in enumerate( variables , 1 ) :

        # y_{k,0}, k > 0, is always true
        p.add_clause( ( y( k , 0 ) , ) )

        # y_{k,t} with 0 < k < t is always false
        for t in range(k+1, pivot+1 ):
            p.add_clause( ( -y( k , t ) , ) )

        for t in range(1,pivot+1):
            p.add_clause( ( -y( k , t ) , y( k - 1 , t ) , y( k - 1 , t - 1 ) ) )
            p.add_clause( ( -y( k , t ) , y( k - 1 , t ) , var ) )

    # we want y_{nvars,pivot} to be true
    p.add_clause( ( y( nvars , pivot ) , ) )

    # return last y + 1
    return y( nvars , pivot ) + 1

def sat_skip ( solution ) :

    lastx = x(R-1,M-1,N-1)

    for i in range( 1 , lastx + 1 ) :

        if solution[i] :
            yield -i

        else :
            yield i


def puzzle ( pivot ) :

    p = SAT( solver = 'cryptominisat' )

    # add increasing constraint
    for k in range( R ) :
        for i1 in range( M ) :
            for j1 in range( N ) :
                for i2 in range( i1 , M ) :
                    for j2 in range( j1 , -1 , -1 ) :
                        if i1 == i2 and j1 == j2 : continue
                        p.add_clause( ( -x(k , i1 , j1) , -x(k , i2 , j2) ) )
                for i2 in range( i1 , -1 , -1 ) :
                    for j2 in range( j1 , N ) :
                        if i1 == i2 and j1 == j2 : continue
                        p.add_clause( ( -x(k , i1 , j1) , -x(k , i2 , j2) ) )

    # add killing constraint
    for k1 in range( R ) :
        for i1 in range( M ) :
            for j1 in range( N ) :
                for k2 in range( k1 + 1 , R ) :
                    p.add_clause( ( -x(k1 , i1 , j1) , -x(k2, i1 , j1 ) ) )
                    for i2 in range( i1 + 1 , M ) :
                        p.add_clause( ( -x(k1, i1, j1) , -x(k2, i2, j1 ) ) )
                    for j2 in range( j1 + 1 , N ) :
                        p.add_clause( ( -x(k1, i1, j1) , -x(k2, i1, j2 ) ) )

    lastx = x(R-1,M-1,N-1)

    max_vars_sat ( p , lastx , lastx+1 , range(1,lastx+1) , pivot )

    return p

def solve_puzzle( pivot ) :
    return puzzle( pivot )( )

LB = min( M , N ) # lower bound
UB = int( M * N ) # upper bound

lb = LB
ub = lb # because of exponential search

print( 'Start exponential search')
while True :

    print( 'ES target {} best {}'.format(lb , best ))

    sol = solve_puzzle(ub)

    if sol :
        output_solution(ub,sol)
        best = ub
        best_sol = sol
        lb = ub + 1
        ub *= 2
    else :
        output_unsat(ub)
        ub = ub - 1
        break

print( 'Start binary search')
while lb <= ub :

    pivot = ( lb + ub ) // 2

    print( 'BS lb {} ub {} pivot {} best {}'.format(lb , ub , pivot , best ))

    sol = solve_puzzle(pivot)

    if sol :
        output_solution(pivot,sol)
        best = pivot
        best_sol = sol
        lb = pivot + 1
    else :
        output_unsat(pivot)
        ub = pivot - 1

print( 'Found one optimal solution with score {}'.format( best ) )

if args.enumerate :

    print( 'Start optimal solutions enumeration' )

    opt = best
    solutions = [ best_sol ]

    while True :

        p = puzzle( opt )

        for solution in solutions :

            p.add_clause( tuple( sat_skip( solution ) ) )

        sol = p( )

        if not sol :

            print( 'No more optimal solutions' )
            break

        output_solution(opt, sol)
        solutions.append( sol )

    print( 'Found {} optimal solution with score {}.'.format( len( solutions ) , opt) )
