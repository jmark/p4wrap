# ifndef P4WRAP_H
# define P4WRAP_H

# if defined(P4_TO_P8)
# include <p8est.h>
# include <p8est_ghost.h>
# include <p8est_extended.h>
# include <p8est_geometry.h>
# include <p8est_connectivity.h>

# define p4wrap_types_mod p8wrap_types_mod
# define p4est_interfaces_mod p8est_interfaces_mod
# define p4wrap_interfaces_mod p8wrap_interfaces_mod

# define p4wrap_mesh_t p8wrap_mesh_t
# define p4wrap_quad_t p8wrap_quad_t
# define p4wrap_side_t p8wrap_side_t

# define p4wrap_gather_proc_offsets p8wrap_gather_proc_offsets
# define p4wrap_gather_mirror_proc_offsets p8wrap_gather_mirror_proc_offsets
# define p4wrap_gather_mirror_proc_mirrors p8wrap_gather_mirror_proc_mirrors

# define p4wrap_quadrant_child_id p8wrap_quadrant_child_id

# define p4wrap_gather_mirroridx p8wrap_gather_mirroridx
# define p4wrap_gather_mirrorptr p8wrap_gather_mirrorptr
# define p4wrap_gather_connection_info p8wrap_gather_connection_info

# define p4wrap_copy_gfq p8wrap_copy_gfq

# define p4wrap_refine p8wrap_refine
# define p4wrap_balance p8wrap_balance
# define p4wrap_coarsen p8wrap_coarsen
# define p4wrap_partition p8wrap_partition

# define P4WRAP_FLAG_NONE P8WRAP_FLAG_NONE
# define P4WRAP_FLAG_REFINE P8WRAP_FLAG_REFINE
# define P4WRAP_FLAG_COARSEN P8WRAP_FLAG_COARSEN
# define P4WRAP_FLAG_ISGHOST P8WRAP_FLAG_ISGHOST

# define p4wrap_get_num_ghosts p8wrap_get_num_ghosts
# define p4wrap_get_num_mirrors p8wrap_get_num_mirrors
# define p4wrap_get_num_quads p8wrap_get_num_quads
# define p4wrap_get_num_global p8wrap_get_num_global
# define p4wrap_get_idx_global p8wrap_get_idx_global
# define p4wrap_get_num_sides p8wrap_get_num_sides
# define p4wrap_get_num_edges p8wrap_get_num_edges
# define p4wrap_get_num_verts p8wrap_get_num_verts
# define p4wrap_get_num_verts_hanging p8wrap_get_num_verts_hanging

# define PP_N_DIMS 3
# define P4EST_CORNERS 8

# else

# include <p4est.h>
# include <p4est_ghost.h>
# include <p4est_extended.h>
# include <p4est_geometry.h>
# include <p4est_connectivity.h>

# define PP_N_DIMS 2
# define P4EST_CORNERS 4
# endif

# define NOR 0
# define SOU 1
# define WES 2
# define EAS 3

# define NWF 0
# define SWF 1
# define NEF 2
# define SEF 3

# define NWB 4
# define SWB 5
# define NEB 6
# define SEB 7

# define p4est_iter_vert_t p4est_iter_corner_t

extern int rank;

enum
{
    P4WRAP_FLAG_NONE    = 0,
    P4WRAP_FLAG_REFINE  = 1,
    P4WRAP_FLAG_COARSEN = 2,
    P4WRAP_FLAG_ISGHOST = 4,
};

typedef struct p4wrap_quad {
    /* Do not change order! */
    int32_t     morton[PP_N_DIMS];
    int8_t      level;
 
    size_t      cellid;
    int32_t     treeid;
    int32_t     quadid;

    size_t      isides[P4EST_FACES];
    size_t      iverts[P4EST_CORNERS];

    int16_t     flags;

    int32_t     local_num;
} p4wrap_quad_t;

typedef struct p4wrap_side_struct {
    int8_t tree_boundary;
    int8_t faceid[2];
    int8_t nquads[2];
    size_t iquads[2][P4EST_HALF];
} p4wrap_side_t;

typedef struct p4wrap_vert_struct {
    int8_t tree_boundary;
    int8_t nquads; // number of touching quadrants
    size_t iquads[P4EST_CORNERS]; // indices to touching quadrants
} p4wrap_vert_t;

typedef struct p4wrap_edge_struct {
    int8_t tree_boundary;
    int8_t nquads[4];
    size_t iquads[4][2];
} p4wrap_edge_t;

# endif
