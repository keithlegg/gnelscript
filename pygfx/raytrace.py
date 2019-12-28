"""
    ORIGINAL CODE IS FROM :
       https://gist.github.com/rossant/6046463

    This is entirely seperate from the other code. 
    I have plans to integrate it into the other code slowly
    It is fun to play with and learn from.

"""


import math 

import numpy as np


#required for framebuffer 
from gnelscript.pygfx.raster_ops import *

#not required, but FUN! 
#from gnelscript.pygfx.point_ops import *
#from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.obj3d import  *

export = object3d() 


class raytracer(object):

    def __init__(self):
        # Default light and material parameters.
        self.ambient     = .1 #.05
        self.diffuse_c   = 1.
        self.specular_c  = 1.
        self.specular_k  = 50
        self.color_light = np.ones(3)

        self.w = 100
        self.h = 100

    def normalize(self, x):
        x /= np.linalg.norm(x)
        return x

    def intersect_plane(self, O, D, P, N):
        # Return the distance from O to the intersection of the ray (O, D) with the 
        # plane (P, N), or +inf if there is no intersection.
        # O and P are 3D points, D and N (normal) are normalized vectors.
        denom = np.dot(D, N)
        if np.abs(denom) < 1e-6:
            return np.inf
        d = np.dot(P - O, N) / denom
        if d < 0:
            return np.inf
        return d

    def intersect_sphere(self, O, D, S, R):
        # Return the distance from O to the intersection of the ray (O, D) with the 
        # sphere (S, R), or +inf if there is no intersection.
        # O and S are 3D points, D (direction) is a normalized vector, R is a scalar.

        a = np.dot(D, D)
      
        OS = O - S
        b = 2 * np.dot(D, OS)
        c = np.dot(OS, OS) - R * R
        disc = b * b - 4 * a * c
        if disc > 0:
            
            distSqrt = np.sqrt(disc) #original 

            # distSqrt = math.sin(np.sqrt(disc)*5) # makes a neat wavy effect 
            # distSqrt = math.cos(np.sqrt(disc)*5) # makes a neat wavy effect 

            q = (-b - distSqrt) / 2.0 if b < 0 else (-b + distSqrt) / 2.0
            t0 = q / a
            t1 = c / q
            t0, t1 = min(t0, t1), max(t0, t1)
            if t1 >= 0:
                return t1 if t0 < 0 else t0
        return np.inf

    def intersect(self, O, D, obj):
        if obj['type'] == 'plane':
            return self.intersect_plane(O, D, obj['position'], obj['normal'])
        elif obj['type'] == 'sphere':
            return self.intersect_sphere(O, D, obj['position'], obj['radius'])

    def get_normal(self, obj, M):
        # Find normal.
        if obj['type'] == 'sphere':
            N = self.normalize(M - obj['position'])
        elif obj['type'] == 'plane':
            N = obj['normal']
        return N
        
    def get_color(self, obj, M):
        color = obj['color']
        if not hasattr(color, '__len__'):
            color = color(M)
        return color

    def get_transparency(self, obj, M):
        """ keith attempting to do a thing """
        transparency = obj['transparency']
        return transparency

    def trace_ray(self, scene, rayO, rayD, L, O):
        # Find first point of intersection with the scene.
        t = np.inf
        for i, obj in enumerate(scene):
            t_obj = self.intersect(rayO, rayD, obj)
            if t_obj < t:
                t, obj_idx = t_obj, i
        # Return None if the ray does not intersect any object.
        if t == np.inf:
            return
        # Find the object.
        obj = scene[obj_idx]
        # Find the point of intersection on the object.
        M = rayO + rayD * t
        
        ##########################
        #keith is having fun here - beware  
        #export.one_vec_to_obj( r3, pos=None, arrowhead=False):        
        #export.one_vec_to_obj( M )
        ##########################


        # Find properties of the object.
        N = self.get_normal(obj, M)

        color = self.get_color(obj, M)
        
        ###########
        # keith attempting to add transparency  
        transparency = self.get_transparency(obj, M)
        if transparency:
            #    print(obj['type'] , 'has transparency ', transparency)
            export.one_vec_to_obj( M )
            if N.any():
                export.one_vec_to_obj( (N[0],N[1],N[2]), M )
        ###########

        toL = self.normalize(L - M)
        toO = self.normalize(O - M)

        # Shadow: find if the point is shadowed or not.
        l = [self.intersect(M + N * .0001, toL, obj_sh) 
                for k, obj_sh in enumerate(scene) if k != obj_idx]
        if l and min(l) < np.inf:
            return
        
        # Start computing the color.
        col_ray = self.ambient
        
        # Lambert shading (diffuse).
        col_ray += obj.get('diffuse_c', self.diffuse_c) * max(np.dot(N, toL), 0) * color
        
        # Blinn-Phong shading (specular).
        col_ray += obj.get('specular_c', self.specular_c) * max(np.dot(N, self.normalize(toL + toO)), 0) ** self.specular_k * self.color_light
        
        return obj, M, N, col_ray

    def add_sphere(self, position, radius, color):
        return dict(type='sphere', position=np.array(position), 
            radius=np.array(radius), color=np.array(color), reflection=0, transparency=.2)
        
    def add_plane(self, position, normal, cp1, cp2):
        return dict(type='plane', position=np.array(position), 
            normal=np.array(normal),
            color=lambda M: (cp1 
                if (int(M[0] * 2) % 2) == (int(M[2] * 2) % 2) else cp2),
            diffuse_c=.75, specular_c=.5, reflection=1, transparency=None)
    
    
    def main(self):
        # List of objects.
        color_plane0 = 1. * np.ones(3)
        color_plane1 = 0. * np.ones(3)

        ## scene = [self.add_sphere([.75, .1, 1.], .6, [0., 0., 1.]),
        ##          self.add_sphere([-.75, .1, 2.25], .6, [.5, .223, .5]),
        ##          self.add_sphere([-2.75, .1, 3.5], .6, [1., .572, .184]),
        ##          self.add_plane([0., -.5, 0.], [0., 1., 0.], color_plane0, color_plane1),
        ##     ]


        scene = [ self.add_sphere([ 0 , 3 , 3 ], 1 , [1. , .572 , .184] ),
                  self.add_sphere([ -4 , 3 , 10 ], 1, [1. , .572 , .184] ),
                  self.add_plane( [0., -5, 0.], [0., 1., 0.], color_plane0, color_plane1),
            ]

        # Light position and color.
        L = np.array([5., 5., -10.])

        depth_max = 5      # Maximum number of light reflections.
        col = np.zeros(3)  # Current color.
        
        O = np.array([0., 0.95, -1.])  # Camera.
        Q = np.array([0., 0., 0.])  # Camera pointing to.

        img = np.zeros((self.h, self.w, 3))

        r = float(self.w) / self.h
        # Screen coordinates: x0, y0, x1, y1.
        S = (-1., -1. / r + .25, 1., 1. / r + .25)

        # Loop through all pixels.
        for i, x in enumerate(np.linspace(S[0], S[2], self.w)):
           
            if i % 10 == 0:
                print( i / float(self.w) * 100, "%" )

            for j, y in enumerate(np.linspace(S[1], S[3], self.h)):
                col[:] = 0
                Q[:2] = (x, y)
                D = self.normalize(Q - O)
                depth = 0
                rayO, rayD = O, D
                reflection = 1.
                # Loop through initial and secondary rays.
                while depth < depth_max:
                    traced = self.trace_ray(scene, rayO, rayD, L, O)
                    if not traced:
                        break
                    obj, M, N, col_ray = traced
                    # Reflection: create a new ray.
                    rayO, rayD = M + N * .0001, self.normalize(rayD - 2 * np.dot(rayD, N) * N)
                    depth += 1
                    col += reflection * col_ray
                    reflection *= obj.get('reflection', 1.)
                img[self.h - j - 1, i, :] = np.clip(col, 0, 1)


        export.save("dump_geom_raytrace.obj") 
        return img 

    def save_image(self, pixdata):
        fb = PixelOp()   
        #fb = RasterObj()
        fb.create_buffer(self.w, self.h)

        fb.insert_numpy( pixdata )
        fb.save('raytraced.png')






