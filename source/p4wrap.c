# include <stdio.h>
# include <stdint.h>

# include "p4wrap.h"

int
p4wrap_quadrant_child_id (const int8_t level, const int32_t x, const int32_t y
# if defined(P4_TO_P8)
    ,const int32_t z
# endif
) {
  int id = 0;

  if (level == 0) {
    return 0;
  }

  id |= ((x & P4EST_QUADRANT_LEN (level)) ? 0x01 : 0);
  id |= ((y & P4EST_QUADRANT_LEN (level)) ? 0x02 : 0);
# if defined(P4_TO_P8)
  id |= ((z & P4EST_QUADRANT_LEN (level)) ? 0x04 : 0);
# endif

  return id;
}

static void _gather_number_sides_cb(p4est_iter_face_info_t *info, void *user_data)
{
    (*(size_t*) user_data)++;
}
 
static size_t _get_local_num_sides(p4est_t* p4est, p4est_ghost_t *layer)
{
    size_t seqid = 0;
    p4est_iterate_ext(p4est, layer, &seqid, NULL, _gather_number_sides_cb,
# if defined(P4_TO_P8) 
        NULL, 
# endif
    NULL, 0);
    return seqid;
}

# if defined(P4_TO_P8) 
// static void _gather_number_edges_cb(p8est_iter_edge_info_t *info, void *user_data)
// {
//     (*(size_t*) user_data)++;
// }
//  
// static size_t _get_local_num_edges(p8est_t* p4est, p8est_ghost_t *layer)
// {
//     size_t seqid = 0;
//     p8est_iterate_ext(p4est,layer,&seqid,NULL,NULL,_gather_number_edges_cb,NULL,0);
//     return seqid;
// }
# endif

static void _gather_number_verts_hanging_face_cb(p4est_iter_face_info_t *info, void *user_data)
{
    for (int i = 0; i < info->sides.elem_count; i++) {
        const p4est_iter_face_side_t *src = p4est_iter_fside_array_index_int(&(info->sides), i);
        if (src->is_hanging) {
            (*(size_t*) user_data)++;
            return;
        }
    }
}

# if defined(P4_TO_P8) 
static void _gather_number_verts_hanging_edge_cb(p8est_iter_edge_info_t *info, void *user_data)
{
    for (int i = 0; i < info->sides.elem_count; i++) {
        const p8est_iter_edge_side_t *src = p8est_iter_eside_array_index_int(&(info->sides),i);
        if (src->is_hanging) {
            (*(size_t*) user_data)++;
            return;
        }
    }
}
# endif

static size_t _get_local_num_verts_hanging(p4est_t* p4est, p4est_ghost_t *layer)
{
    size_t seqid = 0;
    p4est_iterate_ext(p4est,layer,&seqid,NULL,_gather_number_verts_hanging_face_cb,
# if defined(P4_TO_P8) 
        _gather_number_verts_hanging_edge_cb,
# endif
    NULL, 0);
    return seqid;
}

static void _gather_number_verts_conforming_cb(p4est_iter_corner_info_t *info, void *user_data)
{
    (*(size_t*) user_data)++;
}
 
static size_t _get_local_num_verts_conforming(p4est_t* p4est, p4est_ghost_t *layer)
{
    size_t seqid = 0;
    p4est_iterate_ext(p4est, layer, &seqid, NULL, NULL,
# if defined(P4_TO_P8) 
        NULL, 
# endif
    _gather_number_verts_conforming_cb, 0);
    return seqid;
}

size_t p4wrap_get_num_ghosts(p4est_ghost_t *layer)
{
    return layer->ghosts.elem_count;
}

size_t p4wrap_get_num_mirrors(p4est_ghost_t *layer)
{
    return layer->mirrors.elem_count;
}

size_t p4wrap_get_num_sides(p4est_t* p4est, p4est_ghost_t *layer)
{
    return _get_local_num_sides(p4est,layer);
}

size_t p4wrap_get_num_verts(p4est_t* p4est, p4est_ghost_t *layer)
{
    return _get_local_num_verts_conforming(p4est,layer);
}

size_t p4wrap_get_num_verts_hanging(p4est_t* p4est, p4est_ghost_t *layer)
{
    return _get_local_num_verts_hanging(p4est,layer);
}

// # if defined(P4_TO_P8) 
// size_t p4wrap_get_num_edges(p8est_t* p4est, p8est_ghost_t *layer)
// {
//     return _get_local_num_edges(p4est,layer);
// }
// # endif

size_t p4wrap_get_num_quads(p4est_t* p4est)
{
    return p4est->local_num_quadrants;
}

size_t p4wrap_get_num_global(p4est_t* p4est)
{
    return p4est->global_num_quadrants;
}

size_t p4wrap_get_idx_global(p4est_t* p4est)
{
    return p4est->global_first_quadrant[p4est->mpirank];
}

void
p4wrap_gather_proc_offsets(const p4est_ghost_t *layer, const int len, size_t *proc_offsets)
{    
    // Use memcopy in the future.
    for (size_t i = 0; i < len; i++) {
        proc_offsets[i] = layer->proc_offsets[i];
    }
}

void
p4wrap_gather_mirror_proc_offsets(const p4est_ghost_t *layer, const int len, size_t *mirror_proc_offsets)
{    
    // Use memcopy in the future.
    for (size_t i = 0; i < len; i++) {
        mirror_proc_offsets[i] = layer->mirror_proc_offsets[i];
    }
}

void
p4wrap_gather_mirror_proc_mirrors(const p4est_ghost_t *layer, const int len, size_t *mirror_proc_mirrors)
{    
    // Use memcopy in the future.
    for (size_t i = 0; i < len; i++) {
        mirror_proc_mirrors[i] = layer->mirror_proc_mirrors[i];
    }
}

void
p4wrap_gather_mirroridx(p4est_ghost_t *layer, size_t *mirroridx)
{    
    for (size_t q = 0; q < layer->mirrors.elem_count; q++) {
        const p4est_quadrant_t *quad = (p4est_quadrant_t*) sc_array_index(&(layer->mirrors), q);
        mirroridx[q] = quad->p.piggy3.local_num;
    }
}

void
p4wrap_gather_mirrorptr(p4est_ghost_t *layer, char *mirrorptr, size_t size)
{    
    for (size_t q = 1; q < layer->mirrors.elem_count; q++) {
        mirrorptr[q] = mirrorptr[0] + q*size;
    }
}

typedef struct cb_user_data_t {
    p4wrap_quad_t *quads;
    p4wrap_side_t *sides;
    p4wrap_vert_t *verts;
    p4wrap_edge_t *edges;

    size_t sides_seqid;
    size_t verts_seqid;
    size_t edges_seqid;
} cb_user_data_t;

# if defined(P4_TO_P8)
static const int8_t map_to_opposite_vertid[P4EST_CORNERS] = {
    [NWF] = SEB,
    [SWF] = NEB,
    [NEF] = SWB,
    [SEF] = NWB,

    [NWB] = SEF,
    [SWB] = NEF,
    [NEB] = SWF,
    [SEB] = NWF,
};
# else
static const int8_t map_to_opposite_vertid[P4EST_CORNERS] = {
    [NWF] = SEF,
    [SWF] = NEF,
    [NEF] = SWF,
    [SEF] = NWF,
};
# endif

static void _gather_face_connection_info_cb(p4est_iter_face_info_t *info, void *user_data)
{
    cb_user_data_t *udata = (cb_user_data_t*) user_data;
    p4wrap_side_t *side = &(udata->sides[udata->sides_seqid]);

    side->tree_boundary = info->tree_boundary;

    // defaults
    for (int i = 0; i < 2; i++) {
        side->nquads[i] =  0;
        side->faceid[i] = -1;

        for (int8_t j = 0; j < P4EST_HALF; j++)
            side->iquads[i][j] = -1;
    }

    // fill side connection info
    for (int i = 0; i < info->sides.elem_count; i++) {
        const p4est_iter_face_side_t *src = p4est_iter_fside_array_index_int(&(info->sides), i);

        side->faceid[i] = src->face;
        side->nquads[i] = src->is_hanging ? P4EST_HALF : 1;

        if (src->is_hanging) {
            // udata->verts_seqid++;
            for (int8_t j = 0; j < P4EST_HALF; j++) {
                const p4est_quadrant_t *quad = src->is.hanging.quad[j];
                const size_t quadid = quad->p.piggy3.local_num;

                side->iquads[i][j] = quadid;
                udata->quads[quadid].isides[src->face] = udata->sides_seqid;
            }
        } else {
            const p4est_quadrant_t *quad = src->is.full.quad;
            const size_t quadid = quad->p.piggy3.local_num;

            side->iquads[i][0] = quadid;
            udata->quads[quadid].isides[src->face] = udata->sides_seqid;
        }
    }

    udata->sides_seqid++;
}

static void _gather_vert_connection_info_cb(p4est_iter_corner_info_t *info, void *user_data)
{
    cb_user_data_t *udata = (cb_user_data_t*) user_data;
    p4wrap_vert_t *vert = &(udata->verts[udata->verts_seqid]);

    vert->tree_boundary = info->tree_boundary;

    vert->nquads = info->sides.elem_count;

    // defaults
    for (int i = 0; i < P4EST_CORNERS; i++) {
        vert->iquads[i] = -1;
    }

    for (int i = 0; i < info->sides.elem_count; i++) {
        const p4est_iter_corner_side_t *corner  = p4est_iter_cside_array_index_int(&(info->sides),i);
        const size_t iquad = corner->quad->p.piggy3.local_num;
        vert->iquads[map_to_opposite_vertid[corner->corner]] = iquad;
        udata->quads[iquad].iverts[corner->corner] = udata->verts_seqid;
    }

    udata->verts_seqid++;
}

// # if defined(P4_TO_P8) 
// static void _gather_edge_connection_info_cb(p8est_iter_edge_info_t *info, void *user_data)
// {
//     cb_user_data_t *udata = (cb_user_data_t*) user_data;
//     p4wrap_edge_t *edge = &(udata->edges[udata->edges_seqid]);
// 
//     edge->tree_boundary = info->tree_boundary;
// 
//     edge->nquads = info->sides.elem_count;
// 
//     // defaults
//     for (int i = 0; i < 4; i++) {
//         edge->nquads[i] =  0;
//         edge->edgeid[i] = -1;
// 
//         for (int8_t j = 0; j < 2; j++)
//             edge->iquads[i][j] = -1;
//     }
// 
//     // fill side connection info
//     for (int i = 0; i < info->sides.elem_count; i++) {
//         const p8est_iter_edge_side_t *src = p8est_iter_eside_array_index_int(&(info->sides),i);
// 
//         edge->nquads[i] = src->is_hanging ? 2 : 1;
// 
//         if (src->is_hanging) {
//             // udata->verts_seqid++;
//             for (int8_t j = 0; j < 2; j++) {
//                 const p8est_quadrant_t *quad = src->is.hanging.quad[j];
// 
//                 const size_t iquad = quad->p.piggy3.local_num;
//                 edge->iquads[map_to_opposite_vertid[corner->corner]][j] = iquad;
//                 udata->quads[iquad].iverts[corner->corner] = udata->verts_seqid;
//             }
//         } else {
//             const p8est_quadrant_t *quad = src->is.full.quad;
//             const size_t iquad = quad->p.piggy3.local_num;
// 
//             edge->iquads[i][0] = quadid;
//             udata->quads[quadid].isides[src->face] = udata->sides_seqid;
//         }
//     }
// 
//     udata->edges_seqid++;
// }
// # endif

void
p4wrap_gather_connection_info(p4est_t *p4est, p4est_ghost_t* layer, 
    p4wrap_quad_t *quads, p4wrap_side_t *sides, p4wrap_vert_t *verts)
{
    size_t seqid = 0;

    // update inner cells
    for (size_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; treeid++) {
        p4est_tree_t *tree = p4est_tree_array_index(p4est->trees, treeid);

        for (size_t q = 0; q < tree->quadrants.elem_count; q++) {
            p4est_quadrant_t *_quad  = sc_array_index(&(tree->quadrants), q);
            p4wrap_quad_t *quad = &quads[seqid];

            quad->morton[0] = _quad->x;
            quad->morton[1] = _quad->y;
# if defined(P4_TO_P8)
            quad->morton[2] = _quad->z;
# endif

            quad->level     = _quad->level;

            quad->quadid    = seqid;
            quad->quadid    = q;
            quad->treeid    = treeid;
        
            quad->flags     = 0;

            for (int8_t k = 0; k < P4EST_FACES; k++) {
                quad->isides[k] = -1;
            }

            for (int8_t k = 0; k < P4EST_CORNERS; k++) {
                quad->iverts[k] = -1;
            }

            _quad->p.piggy3.local_num = seqid++;
            _quad->p.piggy3.which_tree = treeid;
        }
    }

    // update ghost cells
    for (size_t q = 0; q < layer->ghosts.elem_count; q++) {
        p4est_quadrant_t *_quad = sc_array_index(&(layer->ghosts), q);
        p4wrap_quad_t *quad = &quads[seqid];

        quad->morton[0] = _quad->x;
        quad->morton[1] = _quad->y;
# if defined(P4_TO_P8)
        quad->morton[2] = _quad->z;
# endif

        quad->level     = _quad->level;

        quad->quadid    = seqid;
        quad->quadid    = q;
        quad->treeid    = -1;

        quad->flags     = 0;
        quad->flags     |= P4WRAP_FLAG_ISGHOST;

        quad->local_num = _quad->p.piggy3.local_num;

        for (int8_t k = 0; k < P4EST_FACES; k++) {
            quad->isides[k] = -1;
        }

        for (int8_t k = 0; k < P4EST_CORNERS; k++) {
            quad->iverts[k] = -1;
        }

        _quad->p.piggy3.local_num = seqid++;
        _quad->p.piggy3.which_tree = -1;
    }

    cb_user_data_t udata;

    udata.quads = quads;
    udata.sides = sides;
    udata.verts = verts;
    //udata.eddge = edges;

    udata.sides_seqid = 0;
    udata.verts_seqid = 0;
    //udata.edges_seqid = 0;

    p4est_iter_face_t iter_face_cb = NULL;
    p4est_iter_vert_t iter_vert_cb = NULL;
# if defined(P4_TO_P8)
    p8est_iter_edge_t iter_edge_cb = NULL;
# endif

    if (sides) {
        iter_face_cb = _gather_face_connection_info_cb;
    }

    if (verts) {
        iter_vert_cb = _gather_vert_connection_info_cb;
    }

// # if defined(P4_TO_P8)
//     if (edges) {
//         iter_edges_cb = _gather_edge_connection_info_cb;
//     }
// # endif

    p4est_iterate_ext(p4est,layer,&udata,NULL,iter_face_cb,
# if defined(P4_TO_P8)
        iter_edge_cb,
# endif
    iter_vert_cb,0);
}

// ========================================================================= //

static int
_probe_refine_cb(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quad)
{
    const int8_t flag = ((int8_t *) p4est->user_pointer)[quad->p.piggy3.local_num];
    if (flag & P4WRAP_FLAG_REFINE) {
        return 1;
    } else {
        return 0;
    }
}

static int
_probe_coarsen_cb(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quads[])
{
    const int8_t flag = ((int8_t *) p4est->user_pointer)[quads[0]->p.piggy3.local_num];
    if (flag & P4WRAP_FLAG_COARSEN) {
        return 1;
    } else {
        return 0;
    }
}

void
p4wrap_refine(p4est_t* p4est, int8_t *flags)
{
    const int non_recursive = 0;
    const int refine_to_maxlevel = -1;
    const p4est_init_t init_cb = NULL;
    const p4est_replace_t replace_cb = NULL;

    void *save_and_restore = p4est->user_pointer;
    p4est->user_pointer = (void*) flags;
    p4est_refine_ext(p4est, non_recursive, refine_to_maxlevel, _probe_refine_cb, init_cb, replace_cb);
    p4est->user_pointer = save_and_restore;
}

void
p4wrap_coarsen(p4est_t* p4est, int8_t *flags)
{
    const int non_recursive = 0;
    const int callback_orphans = 0;
    const p4est_init_t init_cb = NULL;
    const p4est_replace_t replace_cb = NULL;

    void *save_and_restore = p4est->user_pointer;
    p4est->user_pointer = (void*) flags;
    p4est_coarsen_ext(p4est, non_recursive, callback_orphans, _probe_coarsen_cb, init_cb, replace_cb);
    p4est->user_pointer = save_and_restore;
}

void
p4wrap_balance(p4est_t* p4est)
{
    // P4EST_CONNECT_FACE   = 21,
    // P4EST_CONNECT_CORNER = 22,
    // P4EST_CONNECT_FULL   = P4EST_CONNECT_CORNER

    const p4est_init_t init_cb = NULL;
    const p4est_replace_t replace_cb = NULL;

    p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_cb, replace_cb);
}

void
p4wrap_partition(p4est_t* p4est)
{
    const int allow_for_coarsening = 1;
    const p4est_weight_t weight_cb = NULL;
    p4est_partition_ext(p4est, allow_for_coarsening, weight_cb);
}

void
p4wrap_copy_gfq(p4est_t* p4est, int count, int64_t *gfq)
{
    for (int i = 0; i < count; i++)
        gfq[i] = p4est->global_first_quadrant[i];
}

# if 0
void p4wrap_get_morton_parent(p4est_quadrant_t *q, p4est_topidx_t *morton)
{
    morton[0] = q->x & ~p4est_QUADRANT_LEN (q->level);
    morton[1] = q->y & ~p4est_QUADRANT_LEN (q->level);
#ifdef P4_TO_P8
    morton[2] = q->z & ~p4est_QUADRANT_LEN (q->level);
#endif
}

int8_t p4wrap_quadrant_get_child_id(const p4est_quadrant_t * q)
{
    int id = 0;

    if (q->level == 0)
      return 0;

    id |= ((q->x & p4est_QUADRANT_LEN (q->level)) ? 0x01 : 0);
    id |= ((q->y & p4est_QUADRANT_LEN (q->level)) ? 0x02 : 0);
#ifdef P4_TO_P8
    id |= ((q->z & p4est_QUADRANT_LEN (q->level)) ? 0x04 : 0);
#endif
    return (int8_t) id;
}

enum p4wrap_BOUNDARY
{
    p4wrap_BOUNDARY_NULL = 0,
    p4wrap_BOUNDARY_NOR  = 1,
    p4wrap_BOUNDARY_SOU  = 2,
    p4wrap_BOUNDARY_WES  = 4,
    p4wrap_BOUNDARY_EAS  = 8
};

int8_t p4wrap_quadrant_get_boundary(const p4est_quadrant_t * q)
{
    p4est_qcoord_t mask;
    uint8_t code;

    if (q->level == 0)
      return p4wrap_BOUNDARY_NOR + p4wrap_BOUNDARY_SOU + p4wrap_BOUNDARY_WES + p4wrap_BOUNDARY_EAS;

    mask = 0;
    for (int8_t i = 0; i < q-> level; i++)
        mask |= 1UL << (p4est_QMAXLEVEL-i);

    code = p4wrap_BOUNDARY_NULL;

    if (q->x == 0)
        code = code + p4wrap_BOUNDARY_NOR;
    else if (q->x == mask)
        code = code + p4wrap_BOUNDARY_SOU;

    if (q->y == 0)
        code = code + p4wrap_BOUNDARY_WES;
    else if (q->y == mask)
        code = code + p4wrap_BOUNDARY_EAS;

    //printf("%3d %3d %20d %20d %20d\n", q->level, code, q->x, q->y, mask);

    return code;
}

// ========================================================================= //
// private auxiliary routines

static inline int8_t _quadrant_get_sibling_id(const p4est_quadrant_t * q)
{
    int id = 0;

    if (q->level == 0)
      return 0;

    id |= ((q->x & p4est_QUADRANT_LEN (q->level)) ? 0x01 : 0);
    id |= ((q->y & p4est_QUADRANT_LEN (q->level)) ? 0x02 : 0);
#ifdef P4_TO_P8
    id |= ((q->z & p4est_QUADRANT_LEN (q->level)) ? 0x04 : 0);
#endif
    return (int8_t) id;
}

// ========================================================================= //

void
p4wrap_mesh_init(p4wrap_mesh_t *mesh)
{    
    mesh->minlevel = 0;
    mesh->payload_size = 0;

    mesh->p4est = NULL;
    mesh->conn = NULL;
    mesh->halo = NULL;

    mesh->cells = NULL;
    mesh->sides = NULL;

    mesh->sendbuf = NULL;
    mesh->recvbuf = NULL;
}
# endif

# if 0
// Debug code pasted here for archiving. Gonna be deleted later.

int twoghosts[2][P4EST_HALF] = {0};
int quadid[2][P4EST_HALF] = {0};
p4est_quadrant_t *quads[2][P4EST_HALF] = {0};

int gg[2] = {0};

gg[0] = 0;
gg[1] = 0;

for (int i = 0; i < info->sides.elem_count; i++) {
    p4est_iter_face_side_t *src = p4est_iter_fside_array_index_int(&(info->sides), i);

    if (src->is_hanging) {
        for (int8_t j = 0; j < P4EST_HALF; j++) {
            twoghosts[i][j] = src->is.hanging.is_ghost[j];
            quads[i][j] = src->is.hanging.quad[j];
            quadid[i][j] = src->is.hanging.quadid[j];
            gg[i] += src->is.hanging.is_ghost[j];
        }
    } else {
        twoghosts[i][0] = src->is.full.is_ghost;
        quads[i][0] = src->is.full.quad;
        quadid[i][0] = src->is.full.quadid;
        gg[i] += src->is.full.is_ghost;
    }
}

if (gg[0] > 0 && gg[1] > 0) {
    printf("gather nsides: %d %10ld | %d %d | %d %d %d %d | %d %d %d %d | %10d %10d %10d %10d | %10d %10d %10d %10d | %10p %10p %10p %10p | %10p %10p %10p %10p\n", info->p4est->mpirank, *(size_t*) user_data, gg[0], gg[1],
        twoghosts[0][0], twoghosts[0][1], twoghosts[0][2], twoghosts[0][3],
        twoghosts[1][0], twoghosts[1][1], twoghosts[1][2], twoghosts[1][3],

        quadid[0][0], quadid[0][1], quadid[0][2], quadid[0][3],
        quadid[1][0], quadid[1][1], quadid[1][2], quadid[1][3],

        (void*) quads[0][0], (void*) quads[0][1], (void*) quads[0][2], (void*) quads[0][3],
        (void*) quads[1][0], (void*) quads[1][1], (void*) quads[1][2], (void*) quads[1][3]
    );
    return;
}
# endif
