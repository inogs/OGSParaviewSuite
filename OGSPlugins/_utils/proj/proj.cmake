##############################################
###         PROJ CMAKE FOR PARAVIEW        ###
##############################################

set( PROJ_FOLDER "${CMAKE_CURRENT_LIST_DIR}/src")

if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wswitch -Wshadow -Wunused-parameter -Wmissing-prototypes -Wmissing-declarations -Wformat -Wformat-security")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Wswitch -Wshadow -Wunused-parameter -Wmissing-declarations -Wformat -Wformat-security")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wswitch -Wshadow -Wunused-parameter -Wmissing-prototypes -Wmissing-declarations -Wformat -Wformat-security -Wfloat-conversion -Wc99-extensions -Wc11-extensions")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Wswitch -Wshadow -Wunused-parameter -Wmissing-declarations -Wformat -Wformat-security -Wfloat-conversion")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

set ( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -I${PROJ_FOLDER}") 
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${PROJ_FOLDER}") 

SET(SRC_PROJ
        ${PROJ_FOLDER}/nad_init.c
        ${PROJ_FOLDER}/PJ_aea.c
        ${PROJ_FOLDER}/PJ_aeqd.c
        ${PROJ_FOLDER}/PJ_airy.c
        ${PROJ_FOLDER}/PJ_aitoff.c
        ${PROJ_FOLDER}/PJ_august.c
        ${PROJ_FOLDER}/PJ_axisswap.c
        ${PROJ_FOLDER}/PJ_bacon.c
        ${PROJ_FOLDER}/PJ_bipc.c
        ${PROJ_FOLDER}/PJ_boggs.c
        ${PROJ_FOLDER}/PJ_bonne.c
        ${PROJ_FOLDER}/PJ_calcofi.c
        ${PROJ_FOLDER}/PJ_cart.c
        ${PROJ_FOLDER}/PJ_cass.c
        ${PROJ_FOLDER}/PJ_cc.c
        ${PROJ_FOLDER}/PJ_ccon.c
        ${PROJ_FOLDER}/PJ_cea.c
        ${PROJ_FOLDER}/PJ_chamb.c
        ${PROJ_FOLDER}/PJ_collg.c
        ${PROJ_FOLDER}/PJ_comill.c
        ${PROJ_FOLDER}/PJ_crast.c
        ${PROJ_FOLDER}/PJ_deformation.c
        ${PROJ_FOLDER}/PJ_denoy.c
        ${PROJ_FOLDER}/PJ_eck1.c
        ${PROJ_FOLDER}/PJ_eck2.c
        ${PROJ_FOLDER}/PJ_eck3.c
        ${PROJ_FOLDER}/PJ_eck4.c
        ${PROJ_FOLDER}/PJ_eck5.c
        ${PROJ_FOLDER}/PJ_eqc.c
        ${PROJ_FOLDER}/PJ_eqdc.c
        ${PROJ_FOLDER}/PJ_eqearth.c
        ${PROJ_FOLDER}/PJ_fahey.c
        ${PROJ_FOLDER}/PJ_fouc_s.c
        ${PROJ_FOLDER}/PJ_gall.c
        ${PROJ_FOLDER}/PJ_geoc.c
        ${PROJ_FOLDER}/PJ_geos.c
        ${PROJ_FOLDER}/PJ_gins8.c
        ${PROJ_FOLDER}/PJ_gnom.c
        ${PROJ_FOLDER}/PJ_gn_sinu.c
        ${PROJ_FOLDER}/PJ_goode.c
        ${PROJ_FOLDER}/PJ_gstmerc.c
        ${PROJ_FOLDER}/PJ_hammer.c
        ${PROJ_FOLDER}/PJ_hatano.c
        ${PROJ_FOLDER}/PJ_helmert.c
        ${PROJ_FOLDER}/PJ_hgridshift.c
        ${PROJ_FOLDER}/PJ_horner.c
        ${PROJ_FOLDER}/PJ_igh.c
        ${PROJ_FOLDER}/PJ_isea.c
        ${PROJ_FOLDER}/PJ_imw_p.c
        ${PROJ_FOLDER}/PJ_krovak.c
        ${PROJ_FOLDER}/PJ_labrd.c
        ${PROJ_FOLDER}/PJ_laea.c
        ${PROJ_FOLDER}/PJ_lagrng.c
        ${PROJ_FOLDER}/PJ_larr.c
        ${PROJ_FOLDER}/PJ_lask.c
        ${PROJ_FOLDER}/PJ_latlong.c
        ${PROJ_FOLDER}/PJ_lcca.c
        ${PROJ_FOLDER}/PJ_lcc.c
        ${PROJ_FOLDER}/PJ_loxim.c
        ${PROJ_FOLDER}/PJ_lsat.c
        ${PROJ_FOLDER}/PJ_misrsom.c
        ${PROJ_FOLDER}/PJ_mbt_fps.c
        ${PROJ_FOLDER}/PJ_mbtfpp.c
        ${PROJ_FOLDER}/PJ_mbtfpq.c
        ${PROJ_FOLDER}/PJ_merc.c
        ${PROJ_FOLDER}/PJ_mill.c
        ${PROJ_FOLDER}/PJ_mod_ster.c
        ${PROJ_FOLDER}/PJ_moll.c
        ${PROJ_FOLDER}/PJ_molodensky.c
        ${PROJ_FOLDER}/PJ_natearth.c
        ${PROJ_FOLDER}/PJ_natearth2.c
        ${PROJ_FOLDER}/PJ_nell.c
        ${PROJ_FOLDER}/PJ_nell_h.c
        ${PROJ_FOLDER}/PJ_nocol.c
        ${PROJ_FOLDER}/PJ_nsper.c
        ${PROJ_FOLDER}/PJ_nzmg.c
        ${PROJ_FOLDER}/PJ_ob_tran.c
        ${PROJ_FOLDER}/PJ_ocea.c
        ${PROJ_FOLDER}/PJ_oea.c
        ${PROJ_FOLDER}/PJ_omerc.c
        ${PROJ_FOLDER}/PJ_ortho.c
        ${PROJ_FOLDER}/PJ_patterson.c
        ${PROJ_FOLDER}/PJ_pipeline.c
        ${PROJ_FOLDER}/PJ_poly.c
        ${PROJ_FOLDER}/PJ_putp2.c
        ${PROJ_FOLDER}/PJ_putp3.c
        ${PROJ_FOLDER}/PJ_putp4p.c
        ${PROJ_FOLDER}/PJ_putp5.c
        ${PROJ_FOLDER}/PJ_putp6.c
        ${PROJ_FOLDER}/PJ_qsc.c
        ${PROJ_FOLDER}/PJ_robin.c
        ${PROJ_FOLDER}/PJ_rpoly.c
        ${PROJ_FOLDER}/PJ_sch.c
        ${PROJ_FOLDER}/PJ_sconics.c
        ${PROJ_FOLDER}/PJ_somerc.c
        ${PROJ_FOLDER}/PJ_sterea.c
        ${PROJ_FOLDER}/PJ_stere.c
        ${PROJ_FOLDER}/PJ_sts.c
        ${PROJ_FOLDER}/PJ_tcc.c
        ${PROJ_FOLDER}/PJ_tcea.c
        ${PROJ_FOLDER}/PJ_times.c
        ${PROJ_FOLDER}/PJ_tmerc.c
        ${PROJ_FOLDER}/PJ_tpeqd.c
        ${PROJ_FOLDER}/PJ_unitconvert.c
        ${PROJ_FOLDER}/PJ_urm5.c
        ${PROJ_FOLDER}/PJ_urmfps.c
        ${PROJ_FOLDER}/PJ_vandg.c
        ${PROJ_FOLDER}/PJ_vandg2.c
        ${PROJ_FOLDER}/PJ_vandg4.c
        ${PROJ_FOLDER}/PJ_vgridshift.c
        ${PROJ_FOLDER}/PJ_wag2.c
        ${PROJ_FOLDER}/PJ_wag3.c
        ${PROJ_FOLDER}/PJ_wag7.c
        ${PROJ_FOLDER}/PJ_wink1.c
        ${PROJ_FOLDER}/PJ_wink2.c
        ${PROJ_FOLDER}/proj_etmerc.c
        ${PROJ_FOLDER}/aasincos.c
        ${PROJ_FOLDER}/adjlon.c
        ${PROJ_FOLDER}/bch2bps.c
        ${PROJ_FOLDER}/bchgen.c
        ${PROJ_FOLDER}/biveval.c
        ${PROJ_FOLDER}/dmstor.c
        ${PROJ_FOLDER}/emess.c
        ${PROJ_FOLDER}/geocent.c
        ${PROJ_FOLDER}/geodesic.c
        ${PROJ_FOLDER}/mk_cheby.c
        ${PROJ_FOLDER}/nad_cvt.c
        ${PROJ_FOLDER}/nad_init.c
        ${PROJ_FOLDER}/nad_intr.c
        ${PROJ_FOLDER}/pj_apply_gridshift.c
        ${PROJ_FOLDER}/pj_apply_vgridshift.c
        ${PROJ_FOLDER}/pj_auth.c
        ${PROJ_FOLDER}/pj_ctx.c
        ${PROJ_FOLDER}/pj_fileapi.c
        ${PROJ_FOLDER}/pj_datum_set.c
        ${PROJ_FOLDER}/pj_datums.c
        ${PROJ_FOLDER}/pj_deriv.c
        ${PROJ_FOLDER}/pj_ell_set.c
        ${PROJ_FOLDER}/pj_ellps.c
        ${PROJ_FOLDER}/pj_errno.c
        ${PROJ_FOLDER}/pj_factors.c
        ${PROJ_FOLDER}/pj_fwd.c
        ${PROJ_FOLDER}/pj_gauss.c
        ${PROJ_FOLDER}/pj_gc_reader.c
        ${PROJ_FOLDER}/pj_geocent.c
        ${PROJ_FOLDER}/pj_gridcatalog.c
        ${PROJ_FOLDER}/pj_gridinfo.c
        ${PROJ_FOLDER}/pj_gridlist.c
        ${PROJ_FOLDER}/PJ_healpix.c
        ${PROJ_FOLDER}/pj_init.c
        ${PROJ_FOLDER}/pj_initcache.c
        ${PROJ_FOLDER}/pj_inv.c
        ${PROJ_FOLDER}/pj_list.c
        ${PROJ_FOLDER}/pj_log.c
        ${PROJ_FOLDER}/pj_malloc.c
        ${PROJ_FOLDER}/pj_math.c
        ${PROJ_FOLDER}/pj_mlfn.c
        ${PROJ_FOLDER}/pj_msfn.c
        ${PROJ_FOLDER}/pj_mutex.c
        ${PROJ_FOLDER}/proj_4D_api.c
        ${PROJ_FOLDER}/pj_internal.c
        ${PROJ_FOLDER}/pj_open_lib.c
        ${PROJ_FOLDER}/pj_param.c
        ${PROJ_FOLDER}/pj_phi2.c
        ${PROJ_FOLDER}/pj_pr_list.c
        ${PROJ_FOLDER}/pj_qsfn.c
        ${PROJ_FOLDER}/pj_release.c
        ${PROJ_FOLDER}/pj_strerrno.c
        ${PROJ_FOLDER}/pj_transform.c
        ${PROJ_FOLDER}/pj_tsfn.c
        ${PROJ_FOLDER}/pj_units.c
        ${PROJ_FOLDER}/pj_utils.c
        ${PROJ_FOLDER}/pj_zpoly1.c
        ${PROJ_FOLDER}/proj_mdist.c
        ${PROJ_FOLDER}/proj_rouss.c
        ${PROJ_FOLDER}/rtodms.c
        ${PROJ_FOLDER}/vector1.c
        ${PROJ_FOLDER}/pj_strtod.c
 )
#        ${PROJ_FOLDER}/emess.hpp
#        ${PROJ_FOLDER}/geocent.hpp
#        ${PROJ_FOLDER}/pj_list.hpp
#        ${PROJ_FOLDER}/proj_internal.hpp
#        ${PROJ_FOLDER}/proj_math.hpp
#        ${PROJ_FOLDER}/${CMAKE_CURRENT_BINARY_DIR}/proj_config.hpp


set(HEADERS_PROJ
        ${PROJ_FOLDER}/projects.hpp
        ${PROJ_FOLDER}/proj_api.hpp
        ${PROJ_FOLDER}/proj.hpp
        ${PROJ_FOLDER}/geodesic.hpp
)
