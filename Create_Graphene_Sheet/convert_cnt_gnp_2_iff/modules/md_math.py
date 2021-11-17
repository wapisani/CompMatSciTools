
import math
import numpy as np

cos = math.cos
sin = math.sin


#########################
#### BEST FIT PLANE #####
#########################

# Argument, p: List of points
def fitplane(p):

    c = centerpoint(p)
    columns = len(p)
    a = np.zeros((3, columns))

    for i in range(0, 3):
        for j in range(0, columns):
            a[i][j] = p[j][i] - c[i]

    u, s, v = np.linalg.svd(a)
    normal = u[:, 2]

    return [c, normal]


#########################
######## CENTER #########
#########################

def centerpoint(p):

    x_coords = []
    y_coords = []
    z_coords = []

    for i in p:
        x_coords.append(i[0])
        y_coords.append(i[1])
        z_coords.append(i[2])

    x_center = np.mean(x_coords)
    y_center = np.mean(y_coords)
    z_center = np.mean(z_coords)

    return [x_center, y_center, z_center]


#########################
####### DISTANCE ########
#########################

# calculate the distance between two points in 3d
#

def dist(t=(), s=()):
    square = (t[0] - s[0])**2 + (t[1] - s[1])**2 + (t[2] - s[2])**2
    distance = square**0.5
    return distance


##########################
### MOLECULE FUNCTIONS ###
##########################

class Molecule:

    def __init__(self, a_init={}):
        self.dictionary = a_init
        self.total = len(a_init)
        self.center = self.center()

    #########################
    ### CHANGE ATOM TYPES ###
    #########################

    def change_atom_types(self, masses={}):
        a = self.dictionary

        back_mass = {1.008: 1, 12.011: 2, 12.01: 2, 15.999: 3, 14.007: 4}

        for id in a:
            a[id].typenum = back_mass[masses[a[id].typenum]]

    #########################
    ### GET GREATEST ATOM ###
    ###   TYPE NUMBER     ###
    #########################

    def get_greatest_atom_type(self):
        greatest = 0
        for id in a:
            if a[id].typenum > greatest:
                greatest = a[id].typenum
        return greatest

    #######################
    ###### TRANSLATE ######
    #######################

    def translate(self, t=()):
        a = self.dictionary

        for id in a:
            a[id].x = a[id].x + t[0]
            a[id].y = a[id].y + t[1]
            a[id].z = a[id].z + t[2]

    ##########################
    ### ROTATE ABOUT AXIS ####
    ### PARALLEL TO X-AXIS ###
    ##########################

    def rotate_x(self, theta, c=None):
        a = self.dictionary

        theta = math.radians(theta)

        if c is not None:
            self.translate((-c[0], -c[1], -c[2]))

        r_m = [[1, 0, 0],
               [0, cos(theta), -1 * sin(theta)],
               [0, sin(theta), cos(theta)]]

        for id in a:
            mpd = Matrix().multiply(r_m, [[a[id].x], [a[id].y], [a[id].z]])
            a[id].x = mpd[0][0]
            a[id].y = mpd[1][0]
            a[id].z = mpd[2][0]

        if c is not None:
            self.translate(c)

    ##########################
    ### ROTATE ABOUT AXIS ####
    ### PARALLEL TO Y-AXIS ###
    ##########################

    def rotate_y(self, theta, c=None):
        a = self.dictionary

        theta = math.radians(theta)

        if c is not None:
            self.translate((-c[0], -c[1], -c[2]))

        r_m = [[cos(theta), 0, sin(theta)],
               [0, 1, 0],
               [-1 * sin(theta), 0, cos(theta)]]

        for id in a:
            mpd = Matrix().multiply(r_m, [[a[id].x], [a[id].y], [a[id].z]])
            a[id].x = mpd[0][0]
            a[id].y = mpd[1][0]
            a[id].z = mpd[2][0]

        if c is not None:
            self.translate(c)

    ##########################
    ### ROTATE ABOUT AXIS ####
    ### PARALLEL TO Z-AXIS ###
    ##########################

    def rotate_z(self, theta, c=None):
        a = self.dictionary

        theta = math.radians(theta)

        if c is not None:
            self.translate((-c[0], -c[1], -c[2]))

        r_m = [[cos(theta), -1 * sin(theta), 0],
               [sin(theta), cos(theta), 0],
               [0, 0, 1]]

        for id in a:
            mpd = Matrix().multiply(r_m, [[a[id].x], [a[id].y], [a[id].z]])
            a[id].x = mpd[0][0]
            a[id].y = mpd[1][0]
            a[id].z = mpd[2][0]

        if c is not None:
            self.translate(c)

    ##############################
    ##### ULTIMATE ROTATION ######
    ##############################

    def rotate(self, theta1, theta2, theta3, c=None):
        self.rotate_x(theta1, c)
        self.rotate_y(theta2, c)
        self.rotate_z(theta3, c)

    ########################
    ####### MAX DIST #######
    ########################

    # calculate the max distance between a point and
    # many atom coordinates
    #
    def maxdist(self, p=()):

        # start by saying the max is zero
        max = 0

        a = self.dictionary

        for id in a:

            # get the coordinates for the ith atom
            ai = (a[id].x, a[id].y, a[id].z)
            d = dist(p, ai)                    # calculate the distance
            if dist > max:                     # if the distance is greater than the tentative max...
                max = d                          # ...then make that the new max

        return max

    #########################
    #### FIND THE CENTER ####
    #########################

    # find a point to be the "center" of many atoms
    # this is done by taking the average of all the coordinates
    #
    # argument: dictionary of atom objects using the
    # atom id as a key
    #
    # returns a tuple
    #

    def center(self):

        # initialize running sum
        x_sum = 0.0
        y_sum = 0.0
        z_sum = 0.0

        a = self.dictionary
        total = self.total

        for id in a:
            x_sum = x_sum + a[id].x         # sum all x-coordinates
            y_sum = y_sum + a[id].y         # sum all y-coordinates
            z_sum = z_sum + a[id].z         # sum all z-coordinates

        # return the average
        t = (x_sum / total, y_sum / total, z_sum / total)
        return t

    #########################
    #### FIND THE RADIUS ####
    #########################

    # find a the radius of sphere that encompasses
    # many atoms
    #
    # argument: dictionary of atom objects and a number to
    # increase the radius by in order to "pad" the sphere

    def radius(self, pad):

        c = self.center
        max_dist_from_center = self.maxdist(c)
        return max_dist_from_center + pad


###########################
######## MATRICES #########
###########################

class Matrix:

    def dimensions(self, m=[]):
        columns = len(m[0])
        rows = len(m)
        dims = (rows, columns)
        return dims

    def add(self, m=[], n=[]):
        d = self.dimensions(m)
        m_plus_n = []

        for i in range(0, d[0]):
            m_plus_n.append([])
            for j in range(0, d[1]):
                m_plus_n[i].append(m[i][j] + n[i][j])

        return m_plus_n

    def scalar_multiply(self, s, m=[]):
        d = self.dimensions(m)
        m_times_s = []

        for i in range(0, d[0]):
            m_times_s.append([])
            for j in range(0, d[1]):
                m_times_s[i].append(s * m[i][j])

        return m_times_s

    def multiply(self, m=[], n=[]):
        dm = self.dimensions(m)
        dn = self.dimensions(n)
        m_times_n = []

        for i in range(0, dm[0]):
            m_times_n.append([])
            for j in range(0, dn[1]):
                sum = 0
                for k in range(0, dm[1]):
                    sum = sum + m[i][k] * n[k][j]
                m_times_n[i].append(sum)

        return m_times_n
