#ifndef _FIELD_CL_H
#define _FIELD_CL_H

#ifdef _WITH_OPENCL

#ifdef _WITH_OPENCL
#include "clinfo.h"
#endif

#include "field.h"
#include "clinfo.h"

//! copy back the field to host memory
//! \param[inout] f a field
void CopyfieldtoCPU(field *f);

void CopyfieldtoGPU(field *f);

void update_physnode_cl(field *f, int ie, cl_mem physnode_cl, real *physnode,
			cl_ulong *time,
			cl_uint nwait, cl_event *wait, cl_event *done);
void set_source_CL(field *f, char *sourcename_cl);
void set_buf_to_zero_cl(cl_mem *buf, MacroCell *mcell, field *f,
			cl_uint nwait, cl_event *wait,  cl_event *done);
void dtfield_CL(field *f, real tnow, cl_mem *dtwn_cl,
		cl_uint nwait, cl_event *wait, cl_event *done);
void DGFlux_CL(MacroCell *mcell, field *f, int d, cl_mem *wn_cl,
	       cl_uint nwait, cl_event *wait, cl_event *done);
void DGVolume_CL(MacroCell *mcell, field *f, cl_mem *dtwn_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done);
void ExtractInterface_CL(MacroCell *mcell, field *f, int d, cl_mem wn_cl,
			 cl_uint nwait, cl_event *wait, cl_event *done);
void InsertInterface_CL(MacroCell *mcell, field *f, int d, cl_mem dtwn_cl,
			cl_uint nwait, cl_event *wait, cl_event *done);
void ExtractedDGInterface_CL(MacroFace *mface, field *f,
			     cl_uint nwait, cl_event *wait, cl_event *done);
void ExtractedDGBoundary_CL(MacroFace *mface, field *f, real tnow,
			    cl_uint nwait, cl_event *wait, cl_event *done);
void DGMacroCellInterface_CL(MacroFace *mface, field *f, cl_mem *wn_cl,
			     cl_uint nwait, cl_event *wait, cl_event *done);
void DGBoundary_CL(MacroFace *mface, field *f, cl_mem *wn_cl, real tnow,
		   cl_uint nwait, cl_event *wait, cl_event *done);
void DGSource_CL(MacroCell *mcell, field *f, real tnow, cl_mem *wn_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done);
void DGMass_CL(MacroCell *mcell, field *f,
	       cl_uint nwait, cl_event *wait, cl_event *done);


void dtfield_extract_CL(field *f, real tnow, cl_mem *wn_cl,
			cl_uint nwait, cl_event *wait, cl_event *done);

void show_cl_timing(field *f);
#endif // _WITH_OPENCL

#endif // _FIELD_CL_H
