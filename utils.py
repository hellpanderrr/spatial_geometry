def circumsphere(vertices):
    '''
    @param vertices: an array of simplices
    @returns: radius of circumsphere, coordinates of circumsphere
    @rtype: tuple
    '''

    # based on https://westy31.home.xs4all.nl/Circumsphere/ncircumsphere.htm and
    # https://codereview.stackexchange.com/questions/77593/calculating-the-volume-of-a-tetrahedron

    from scipy.spatial import distance

    # distances
    squared_dists = distance.pdist(vertices, metric='sqeuclidean')
    n = distance.num_obs_y(squared_dists)

    # add ones and a zero to make a Cayley-Menger matrix.
    squared_dists_mat = distance.squareform(squared_dists)
    with_border = np.insert(np.insert(squared_dists_mat, 0, values=1, axis=1), 0, values=1, axis=0)
    np.fill_diagonal(with_border, 0)
    
    #the first diagonal element of its inverse holds -2*r*r
    #and first row / first column hold barycentric coordinates of the sphere
    inv = np.linalg.inv(with_border)
    r = math.sqrt(inv[0][0] / -2)
    barycentric_coodinates = inv[1:, 0]
    return r, bary2cart(vertices, barycentric_coodinates)


def bary2cart(simplex, point):
    simplex, point = np.asarray(simplex), np.asarray(point)
    return np.dot(simplex.T, point)


def check_point_in_sphere(sphere, r, point):
    '''
    @param sphere: coodinates of a shpere
    @type sphere: list
    @param r: raduius of a sphere
    @type r: float
    @param point: coordinates of a point to check
    @type point: list
    @returns: boolean value
    @rtype: bool
    '''
    sphere, point = np.asarray(sphere), np.asarray(point)
    if np.sum((sphere - point)**2) <= r ** 2:
        return True
    else:
        return False
