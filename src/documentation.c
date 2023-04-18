/*
 * Autogenerated file, do not  edit
 * Tue Apr 18 18:13:46 BST 2023
 * -------------------------------
 */


/**
 * @ingroup surfaces
 *
 * @brief allocate an ::nbi_surface_t
 *
 * @param nnmax maximum number of nodes in surface;
 * @param npmax maximum number of patches in surface.
 *
 * @return pointer to newly allocated ::nbi_surface_t
 **/

nbi_surface_t *nbi_surface_alloc(gint nnmax, gint npmax) ;

/**
 * @ingroup surfaces
 *
 * @brief write an ::nbi_surface_t to file
 *
 * @param s surface to write to file
 * @param f output file stream
 *
 * @return 0 on success
 **/

gint nbi_surface_write(nbi_surface_t *s, FILE *f) ;

/**
 * @ingroup surfaces
 *
 * @brief allocate and read an ::nbi_surface_t from file
 *
 * @param f input file stream.
 *
 * @return pointer to newly allocated surface, read from file
 **/

nbi_surface_t *nbi_surface_read(FILE *f) ;

/**
 * @ingroup surfaces
 *
 * @brief estimate centroid of surface patch
 *
 * @param x coordinates of patch nodes;
 * @param xstr stride from first component of one to first component of next;
 * @param w patch quadrature weights;
 * @param wstr quadrature weight stride;
 * @param nx number of nodes;
 * @param c on exit contains coordinates of centroid.
 *
 * @return 0 on success.
 **/

gint nbi_surface_patch_centroid(gdouble *x, gint xstr,
				      gdouble *w, gint wstr,
				      gint nx,
				      gdouble *c) ;

/**
 * @ingroup surfaces
 *
 * @brief estimate radius of circumsphere of surface patch 
 *
 * @param x coordinates of patch nodes;
 * @param xstr stride from first component of one node to first component 
 * of next;
 * @param nx number of patch nodes;
 * @param c centroid of patch (from 
 * nbi_surface_patch_centroid)).
 *
 * @return patch radius
 **/

gdouble nbi_surface_patch_radius(gdouble *x, gint xstr,
					gint nx, gdouble *c) ;

/**
 * @ingroup surfaces
 *
 * @brief identify neighbour nodes of surface patch
 *
 * @param c centroid;
 * @param r radius of neighbourhood;
 * @param x list of candidate neighbour nodes;
 * @param xstr stride from first component of one node to first component 
 * of next;
 * @param nx number of candidate neighbour nodes;
 * @param n0 first entry in \a x to check;
 * @param n1 last entry in \a x to check;
 * @param nbrs on exit contains indices of nodes lying within distance 
 * \a r of \a c; 
 * @param nnbrs on exit contains number of nodes in \a nbrs;
 * @param nnmax maximum number of neighbours to be inserted in \a nbrs.
 *
 * @return 0 on success.
 **/

gint nbi_patch_neighbours(gdouble *c, gdouble r,
				gdouble *x, gint xstr, gint nx,
				gint n0, gint n1,
				gint *nbrs, gint *nnbrs,
				gint nnmax) ;

/**
 * @ingroup surfaces
 *
 * @brief identify neighbours of a surface patch
 *
 * @param s an ::nbi_surface_t;
 * @param p index of patch whose neighbours are to be found;
 * @param r radius of neighbourhood;
 * @param nbrs on exit contains indices of nodes of \a s lying with
 * distance \a r of centroid of patch \a p of \a s;
 * @param nnbrs on exit contains number of entries in \a nbrs;
 * @param nnmax limit on number of entries in \a nbrs.
 *
 * @return 0 on success.
 **/

gint nbi_surface_patch_neighbours(nbi_surface_t *s,
					gint p, gdouble r,
					gint *nbrs, gint *nnbrs,
					gint nnmax) ;

/**
 * @ingroup surfaces
 *
 * @brief Koornwinder interpolation matrix for surface patch
 *
 * @param ns number of interpolation nodes on patch;
 * @param K on exit, pointer to interpolation matrix for patch;
 * @param Nk on exit, maximum order of Koornwinder polynomial in interpolation.
 *
 * @return 0 on success.
 **/

gint nbi_element_interp_matrix(gint ns, gdouble **K, gint *Nk) ;

/**
 * @ingroup surfaces
 *
 * @brief matrix for upsampling patch nodes to quadrature nodes
 *
 * @param ns number of patch nodes (must match an available quadrature rule);
 * @param nu number of upsampled nodes (must match an available quadrature 
 * rule).
 *
 * @return pointer to upsampling matrix.
 **/

gdouble *nbi_patch_upsample_matrix(gint ns, gint nu) ;

/**
 * @ingroup surfaces
 *
 * @brief set basic data (centroid and radius) for each patch of a surface.
 *
 * @param s an ::nbi_surface_t whose data are to be set.
 *
 * @return 0 on success
 **/

gint nbi_surface_set_patch_data(nbi_surface_t *s) ;

/**
 * @ingroup surfaces
 *
 * @brief generate a new ::nbi_matrix_t for use in solvers
 *
 * @param s an ::nbi_surface_t
 *
 * @return pointer to new ::nbi_matrix_t corresponding to \a s.
 **/

nbi_matrix_t *nbi_matrix_new(nbi_surface_t *s) ;

/**
 * @ingroup surfaces
 *
 * @brief read an ::nbi_matrix_t from file
 *
 * @param f input file stream;
 * @param m pointer to an allocated ::nbi_matrix_t.
 *
 * @return 0 on success.
 **/

gint nbi_matrix_read(FILE *f, nbi_matrix_t *m) ;

/**
 * @ingroup matrix
 *
 * @brief write an ::nbi_matrix_t to file
 *
 * @param f output file stream;
 * @param m an ::nbi_matrix_t.
 *
 * @return 0 on success.
 **/

gint nbi_matrix_write(FILE *f, nbi_matrix_t *m) ;

/**
 * @ingroup matrix
 *
 * @brief write numerical data to file.
 *
 * The function writes a header which can be interpreted by
 * ::nbi_data_read and writes \a nd lines of data in the form:
 * 
 *    dat[i*dstr+0] dat[i*dstr+1] ... dat[i*dstr+ne-1] 
 * 
 * with \f$i=0,\ldots,nd-1\f$.
 *
 * @param f output file stream;
 * @param dat data;
 * @param dstr stride between first components of successive data entries;
 * @param ne number of components per data entry;
 * @param nd number of data entries.
 *
 * @return 0 on success.
 **/

gint nbi_data_write(FILE *f, gdouble *dat,
			  gint dstr, gint ne, gint nd) ;

/**
 * @ingroup surfaces
 *
 * @brief read numerical data from a file written with
 * ::nbi_data_write
 *
 * @param f input file stream;
 * @param nd on exit contains number of data entries;
 * @param ne on exit contains number of components per entry.
 *
 * @return pointer to newly allocated block of data.
 **/

gdouble *nbi_data_read(FILE *f, gint *nd, gint *ne) ;

/**
 * @ingroup matrix
 *
 * @brief initialize Fast Multipole Method matrix
 *
 * @param m ::nbi_matrix_t to be initialized;
 * @param problem ::nbi_problem_t (Laplace or Helmholtz);
 * @param shifts WBFMM shift operators;
 * @param order_s order of singular expansion at leaf nodes;
 * @param sstr increment in order of singular expansion at each level;
 * @param order_r order of regular expansion at leaf nodes;
 * @param rstr increment in order of regular expansion at each level;
 * @param depth depth of FMM tree;
 * @param dtree increment in size of FMM box to contain all points;
 * @param shift_bw if TRUE (recommended) use backward shifting in FMM;
 * @param precompute_local if TRUE (not recommended) precompute local
 * interactions in FMM;
 * @param work workspace.
 *
 * @return 0 on success.
 **/

gint nbi_matrix_fmm_init(nbi_matrix_t *m,
			       nbi_problem_t problem,
			       wbfmm_shift_operators_t *shifts,
			       guint *order_s, gint sstr,
			       guint *order_r, gint rstr,
			       gint depth,
			       gdouble dtree,
			       gboolean shift_bw,
			       gboolean precompute_local,
			       gdouble *work) ;

/**
 * @ingroup matrix
 *
 * @brief perform matrix multiplication \f$y := \beta y + \alpha A x\f$
 *
 * @param A ::nbi_matrix_t
 * @param x input vector;
 * @param xstr stride in \a x;
 * @param al weight \f$\alpha\f$;
 * @param y output vector;
 * @param ystr stride in \a y;
 * @param bt weight \f$\beta\f$;
 * @param nthreads number of threads to use (if available) in FMM and near-field
 * corrections;
 * @param work workspace.
 *
 * @return 0 on success.
 **/

gint nbi_matrix_multiply(nbi_matrix_t *A,
			       gdouble *x, gint xstr, gdouble al,
			       gdouble *y, gint ystr, gdouble bt,
			       gint nthreads,
			       gdouble *work) ;

/**
 * @ingroup boundary
 *
 * @brief compute boundary condition on surface
 *
 * @param s ::nbi_surface_t on which to set boundary condition;
 * @param p on exit contains ``surface potential'' boundary condition;
 * @param pstr stride in \a p;
 * @param pn on exit contains ``normal derivative'' boundary condition;
 * @param nstr stride in \a pn;
 * @param b boundary condition to apply.
 *
 * @return 0 on success.
 **/

gint nbi_boundary_condition_set(nbi_surface_t *s,
				      gdouble *p , gint pstr,
				      gdouble *pn, gint nstr,
				      nbi_boundary_condition_t *b) ;
