/*
 * Autogenerated file, do not  edit
 * Mon Dec 11 14:52:30 GMT 2023
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
 * @ingroup matrix
 *
 * @brief generate a new ::nbi_matrix_t for use in solvers
 *
 * @param s an ::nbi_surface_t
 *
 * @return pointer to new ::nbi_matrix_t corresponding to \a s.
 **/

nbi_matrix_t *nbi_matrix_new(nbi_surface_t *s) ;

/**
 * @ingroup matrix
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
 * @ingroup export
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
 * @ingroup export
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

/** 
 * @ingroup matrix
 *
 * Assemble a matrix for the Helmholtz boundary integral problem
 * 
 * @param s surface for problem;
 * @param k wavenumber in Helmholtz equation
 * @param eta separation parameter for adaptive quadrature;
 * @param nqa number of points in adaptive quadrature rule;
 * @param dmax maximum recursion depth for adaptive quadrature;
 * @param tol tolerance for adaptive quadrature;
 * @param N order of singular quadrature rule;
 * @param nu number of upsampled sources per patch;
 * @param nnmax maximum number of neighbours in patch neighbour lists;
 * @param nthreads number of threads.
 * 
 * @return newly allocated and assembled ::nbi_matrix_t
 */

nbi_matrix_t *nbi_matrix_assemble_helmholtz(nbi_surface_t *s,
						  gdouble k,
						  gdouble eta,
						  gint nqa, gint dmax,
						  gdouble tol,
						  gint N, gint nu,
						  gint nnmax, gint nthreads) ;

/** 
 * @ingroup matrix
 *
 * Assemble a matrix for the Laplace boundary integral problem
 * 
 * @param s surface for problem;
 * @param eta separation parameter for adaptive quadrature;
 * @param nqa number of points in adaptive quadrature rule;
 * @param dmax maximum recursion depth for adaptive quadrature;
 * @param tol tolerance for adaptive quadrature;
 * @param N order of singular quadrature rule;
 * @param nu number of upsampled sources per patch;
 * @param nnmax maximum number of neighbours in patch neighbour lists;
 * @param nthreads number of threads.
 * 
 * @return newly allocated and assembled ::nbi_matrix_t
 */

nbi_matrix_t *nbi_matrix_assemble_laplace(nbi_surface_t *s,
						gdouble eta,
						gint nqa, gint dmax,
						gdouble tol,
						gint N, gint nu,
						gint nnmax, gint nthreads) ;

/** 
 * @ingroup export
 * 
 * Generate a triangulation of a surface, by recursive subdivision of
 * patches into plane triangles. Optionally, surface data can also be
 * interpolated onto the triangulation.
 * 
 * @param s an ::nbi_surface_t;
 * @param dmax recursion depth for patch splitting;
 * @param K Koornwinder interpolation matrix;
 * @param Nk order of \a K (see ::nbi_element_interp_matrix);
 * @param x on output contains nodes of triangulation;
 * @param xstr stride of points in \a x (must be three or greater);
 * @param npmax maximum number of points in triangulation;
 * @param tri on exit contains triangles as list of indices into \a x;
 * @param tstr stride of triangles in \a tri (must be three or greater);
 * @param ntmax maximum number of triangles in triangulation;
 * @param f if not NULL, data to be interpolated onto triangulation nodes;s
 * @param fstr stride in \a f;
 * @param fi on exit, contains \a f interpolated onto points \a x;
 * @param fistr stride of data in \a fi;
 * @param np on exit, number of points in triangulation;
 * @param nt on exit, number triangles in triangulation.
 * 
 * @return 0 on success.
 */

gint nbi_mesh_triangulate(nbi_surface_t *s,
				gint dmax,
				gdouble *K, gint Nk,
				gdouble *x, gint xstr,
				gint npmax,
				gint *tri, gint tstr,
				gint ntmax,
				gdouble *f, gint fstr,
				gdouble *fi, gint fistr,
				gint *np, gint *nt) ;

/** 
 * @ingroup export
 *
 * Export a triangulation, from ::nbi_mesh_triangulate to a GMSH .msh file
 * 
 * @param f output file stream;
 * @param view if not NULL, name of view in .msh file; if NULL name of view
 * defaults to "NBI view";
 * @param x points of triangulation;
 * @param xstr stride between points in \a x;
 * @param np number of points in \a x;
 * @param offp offset to be added to the point indices written to file (this
 * can be used to avoid conflicts when multiple meshes are read into GMSH);
 * @param tri indices of triangle points;
 * @param tstr stride between triangles in \a tri;
 * @param nt number of triangles in \a tri;
 * @param offt offset to be added to element indices written to file;
 * @param data if not NULL, node data to be written to file;
 * @param dstr stride in \a data.
 * 
 * @return 0 on success.
 */

gint nbi_mesh_export_gmsh(FILE *f, gchar *view,
				gdouble *x, gint xstr, gint np,
				gint offp,
				gint *tri, gint tstr, gint nt,
				gint offt,
				gdouble *data, gint dstr) ;

/** 
 * @ingroup surfaces
 *
 * Generate a spherical surface
 * 
 * @param s surface to which sphere patches are to be added;
 * @param r radius of sphere;
 * @param nth number of patches in azimuth \f$\theta\f$;
 * @param nph number of patches in elevation \f$\phi\f$;
 * @param nq number of quadrature points per patch.
 * 
 * @return 0 on success.
 */

gint nbi_geometry_sphere(nbi_surface_t *s,
			       gdouble r, gint nth, gint nph,
			       gint nq) ;
/** 
 * @ingroup surfaces
 *
 * Generate an ellipsoidal surface
 * 
 * @param s surface to which patches are to be added;
 * @param a ellipsoid length on \f$x\f$ axis;
 * @param b ellipsoid length on \f$y\f$ axis;
 * @param c ellipsoid length on \f$z\f$ axis;
 * @param nth number of patches in azimuth \f$\theta\f$;
 * @param nph number of patches in elevation \f$\phi\f$;
 * @param nq number of quadrature points per patch.
 * 
 * @return 0 on success.
 */

gint nbi_geometry_ellipsoid(nbi_surface_t *s,
				  gdouble a, gdouble b,
				  gdouble c,
				  gint nth, gint nph,
				  gint nq) ;

/** 
 * @ingroup surfaces
 *
 * Generate a spherical surface by recursive subdivision of an
 * icosahedron
 * 
 * @param s surface to which sphere patches are to be added;
 * @param r radius of sphere;
 * @param nr number of subdivisions;
 * @param nq number of quadrature points per patch.
 * 
 * @return 0 on success.
 */

gint nbi_geometry_sphere_ico(nbi_surface_t *s,
				   gdouble r, gint nr, gint nq) ;

/** 
 * @ingroup surfaces
 *
 * Generate an ellipsoidal surface by recursive subvision of a
 * deformed icosahedron
 * 
 * @param s surface to which patches are to be added;
 * @param a ellipsoid length on \f$x\f$ axis;
 * @param b ellipsoid length on \f$y\f$ axis;
 * @param c ellipsoid length on \f$z\f$ axis;
 * @param nr number of subdivisions;
 * @param nq number of quadrature points per patch.
 * 
 * @return 0 on success.
 */

gint nbi_geometry_ellipsoid_ico(nbi_surface_t *s,
				      gdouble a, gdouble b,
				      gdouble c,
				      gint nr, gint nq) ; 

/** 
 * @ingroup surfaces
 *
 * Generate a regular grid in plane \f$z=0\f$ which can be transformed
 * and used for visualisation. Nodes are interpolated between
 * \f$(s_{\min},t_{\min},0)\f$ and \f$(s_{\max},t_{\max},0)\f$
 * 
 * @param s surface to add patches to;
 * @param smin coordinate of lower left corner of grid;
 * @param smax coordinate of upper right corner of grid;
 * @param ns number of nodes in \f$s\f$.
 * @param tmin coordinate of lower left corner of grid;
 * @param tmax coordinate of upper right corner of grid;
 * @param nt number of nodes in \f$t\f$.
 * @param nq number of quadrature nodes per patch on \a s.
 * 
 * @return 0 on succes.
 */

gint nbi_geometry_grid(nbi_surface_t *s,
			     gdouble smin, gdouble smax,
			     gint ns,
			     gdouble tmin, gdouble tmax,
			     gint nt,
			     gint nq) ;
  
/** 
 * @ingroup matrix
 *
 * Upsample nodal data to upsample nodes, multiplying by real weights
 * 
 * @param m matrix to upsample;
 * @param p complex surface potential to upsample;
 * @param pstr stride in \a p;
 * @param pwt weight for \a p;
 * @param pn complex surface potential gradient to upsample;
 * @param nstr stride in \a pn;
 * @param nwt weight for \a pn.
 * 
 * @return 0 on success.
 */

gint nbi_matrix_upsample_helmholtz(nbi_matrix_t *m,
					 gdouble *p, gint pstr,
					 gdouble pwt,
					 gdouble *pn, gint nstr,
					 gdouble nwt) ;

/** 
 * @ingroup matrix
 * 
 * @brief Evaluate Green's identity for Helmholtz equation on surface. 
 *
 * Evaluate Green's identity for Helmholtz equation on surface.  This
 * is mainly used for checking discretizations and boundary condition
 * evaluation.
 * 
 * @param m matrix for boundary integral problem;
 * @param p complex surface potential to upsample;
 * @param pstr stride in \a p;
 * @param pwt weight for \a p;
 * @param pn complex surface potential gradient to upsample;
 * @param nstr stride in \a pn;
 * @param nwt weight for \a pn.
 * @param f on exit contains surface potential at each node of \a s
 * evaluated using Green's identity;
 * @param fstr stride in \a f;
 * @param nthreads number of threads to use;
 * @param work workspace.
 * 
 * @return 0 on success.
 */

gint nbi_surface_greens_identity_helmholtz(nbi_matrix_t *m,
						 gdouble *p ,
						 gint pstr,
						 gdouble pwt,
						 gdouble *pn,
						 gint nstr,
						 gdouble nwt,
						 gdouble *f,
						 gint fstr,
						 gint nthreads,
						 gdouble *work) ;

/** 
 * @ingroup matrix
 * 
 * @brief Perform matrix multiplication for Helmholtz problem
 * 
 * Perform the matrix multiplication (actually corrected FMM
 * evaluation) corresponding to operation \f$y\to \alpha A x + \beta
 * y\f$ (similar to BLAS implementation), mainly for use in iterative solvers. 
 * 
 * @param A matrix for problem;
 * @param x input vector;
 * @param xstr stride in \a x;
 * @param al \f$\alpha\f$;
 * @param y on output contains \f$y\to \alpha A x + \beta\f$
 * @param ystr stride in \a y;
 * @param bt \f$\beta\f$;
 * @param nthreads number of threads to use;
 * @param work workspace.
 * 
 * @return 0 on success.
 */

gint nbi_matrix_multiply_helmholtz(nbi_matrix_t *A,
					 gdouble *x, gint xstr,
					 gdouble al,
					 gdouble *y, gint ystr,
					 gdouble bt,
					 gint nthreads,
					 gdouble *work) ;

/** 
 * @ingroup matrix
 *
 * @brief Calculate field from surface
 *
 * Calculate field radiated from surface using specified boundary
 * conditions. This is not FMM accelerated. 
 * 
 * @param m problem matrix;
 * @param p complex surface potential;
 * @param pstr stride in \a p;
 * @param pwt weight for \a p;
 * @param pn complex surface potential gradient;
 * @param nstr stride in \a pn;
 * @param nwt weight for \a pn.
 * @param x point for field evaluation;
 * @param f on exit contains \f$\phi(\mathbf{x})\f$;
 * @param nthreads number of threads to use;
 * @param work workspace.
 * 
 * @return 0 on success.
 */

gint nbi_calc_field_helmholtz(nbi_matrix_t *m,
				    gdouble *p ,
				    gint pstr,
				    gdouble pwt,
				    gdouble *pn,
				    gint nstr,
				    gdouble nwt,
				    gdouble *x,
				    gdouble *f,
				    gint nthreads,
				    gdouble *work) ;
/** 
 * @ingroup matrix
 *
 * @brief Calculate field from surface
 *
 * Calculate field radiated from surface using specified boundary
 * conditions packed in a single vector. This is not FMM accelerated.
 * 
 * @param s surface for calculation;
 * @param k wavenumber;
 * @param ps complex surface potential and potential gradientto upsample;
 * @param pstr stride in \a ps;
 * @param al \f$\alpha\f$;
 * @param bt \f$\beta\f$;
 * @param x point for field evaluation;
 * @param p on exit contains \f$\beta p + \alpha \phi(\mathbf{x})\f$, 
 * with \f$\phi\f$ calculated using Green's identity. 
 * 
 * @return 0 on success.
 */

gint nbi_surface_field_helmholtz(nbi_surface_t *s,
				 gdouble k,
				 gdouble *ps,
				 gint pstr,
				 gdouble *al, gdouble *bt,
				 gdouble *x,
				 gdouble *p) ;
/** 
 * @ingroup matrix
 *
 * Upsample nodal data to upsample nodes, multiplying by real weights
 * 
 * @param m matrix to upsample;
 * @param p complex surface potential to upsample;
 * @param pstr stride in \a p;
 * @param pwt weight for \a p;
 * @param pn complex surface potential gradient to upsample;
 * @param nstr stride in \a pn;
 * @param nwt weight for \a pn.
 * 
 * @return 0 on success.
 */

gint nbi_matrix_upsample_laplace(nbi_matrix_t *m,
					 gdouble *p, gint pstr,
					 gdouble pwt,
					 gdouble *pn, gint nstr,
					 gdouble nwt) ;

/** 
 * @ingroup matrix
 * 
 * @brief Evaluate Green's identity for Laplace equation on surface. 
 *
 * Evaluate Green's identity for Laplace equation on surface.  This
 * is mainly used for checking discretizations and boundary condition
 * evaluation.
 * 
 * @param m matrix for boundary integral problem;
 * @param p complex surface potential to upsample;
 * @param pstr stride in \a p;
 * @param pwt weight for \a p;
 * @param pn complex surface potential gradient to upsample;
 * @param nstr stride in \a pn;
 * @param nwt weight for \a pn.
 * @param f on exit contains surface potential at each node of \a s
 * evaluated using Green's identity;
 * @param fstr stride in \a f;
 * @param nthreads number of threads to use;
 * @param work workspace.
 * 
 * @return 0 on success.
 */

gint nbi_surface_greens_identity_laplace(nbi_matrix_t *m,
						 gdouble *p ,
						 gint pstr,
						 gdouble pwt,
						 gdouble *pn,
						 gint nstr,
						 gdouble nwt,
						 gdouble *f,
						 gint fstr,
						 gint nthreads,
						 gdouble *work) ;

/** 
 * @ingroup matrix
 * 
 * @brief Perform matrix multiplication for Laplace problem
 * 
 * Perform the matrix multiplication (actually corrected FMM
 * evaluation) corresponding to operation \f$y\to \alpha A x + \beta
 * y\f$ (similar to BLAS implementation), mainly for use in iterative solvers. 
 * 
 * @param A matrix for problem;
 * @param x input vector;
 * @param xstr stride in \a x;
 * @param al \f$\alpha\f$;
 * @param y on output contains \f$y\to \alpha A x + \beta\f$
 * @param ystr stride in \a y;
 * @param bt \f$\beta\f$;
 * @param nthreads number of threads to use;
 * @param work workspace.
 * 
 * @return 0 on success.
 */

gint nbi_matrix_multiply_laplace(nbi_matrix_t *A, gdouble *x, gint xstr,
				       gdouble al, gdouble *y, gint ystr,
				       gdouble bt, gint nthreads,
				       gdouble *work) ;
/** 
 * @ingroup matrix
 *
 * @brief Calculate field from surface
 *
 * Calculate field radiated from surface using specified boundary
 * conditions packed in a single vector. This is not FMM accelerated.
 * 
 * @param s surface for calculation;
 * @param ps complex surface potential and potential gradientto upsample;
 * @param pstr stride in \a ps;
 * @param al \f$\alpha\f$;
 * @param bt \f$\beta\f$;
 * @param x point for field evaluation;
 * @param p on exit contains \f$\beta p + \alpha \phi(\mathbf{x})\f$, 
 * with \f$\phi\f$ calculated using Green's identity. 
 * 
 * @return 0 on success.
 */

gint nbi_surface_field_laplace(nbi_surface_t *s, gdouble *ps, gint pstr,
			       gdouble al, gdouble bt, gdouble *x,
			       gdouble *p) ;

/** 
 * @ingroup matrix
 * 
 * @brief basic GMRES solver for real problems
 * 
 * This is a basic implementation of a GMRES solver. You may find the
 * PETSC solver more useful, if it is available.
 *
 * @param A problem matrix;
 * @param n matrix size;
 * @param x initial guess for solution;
 * @param xstr stride in \a x;
 * @param b right hand side vector;
 * @param bstr stride in \a b;
 * @param m restart interval;
 * @param max_it maximum number of iterations;
 * @param tol convergence tolerance;
 * @param error on exit, contains error estimate;
 * @param nthreads number of threads;
 * @param work workspace, sized using ::nbi_gmres_workspace_size_real.
 * 
 * @return 0 on success.
 */

gint nbi_gmres_real(nbi_matrix_t *A, gint n,
			  gdouble *x, gint xstr,
			  gdouble *b, gint bstr,
			  gint m, gint max_it,
			  gdouble tol, gdouble *error,
			  gint nthreads,
			  gdouble *work) ;

/** 
 * @ingroup matrix
 * 
 * @brief basic GMRES solver for complex problems
 * 
 * This is a basic implementation of a GMRES solver. You may find the
 * PETSC solver more useful, if it is available.
 * 
 * @param A problem matrix;
 * @param n matrix size;
 * @param x initial guess for solution;
 * @param xstr stride in \a x;
 * @param b right hand side vector;
 * @param bstr stride in \a b;
 * @param m restart interval;
 * @param max_it maximum number of iterations;
 * @param tol convergence tolerance;
 * @param error on exit, contains error estimate;
 * @param nthreads number of threads;
 * @param work workspace, sized using ::nbi_gmres_workspace_size_complex.
 * 
 * @return 0 on success.
 */

gint nbi_gmres_complex(nbi_matrix_t *A, gint n,
			     gdouble *x, gint xstr,
			     gdouble *b, gint bstr,
			     gint m, gint max_it,
			     gdouble tol, gdouble *error,
			     gint nthreads,
			     gdouble *work) ;
