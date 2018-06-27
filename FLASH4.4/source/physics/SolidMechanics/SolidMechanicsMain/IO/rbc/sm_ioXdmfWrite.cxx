#include <Xdmf.h>
// We want the filenames to be based on the iteration
// and padded with zeros
using std::setw;
using std::setfill;

// This works with g77. Different compilers require different
// name mangling.
#define XdmfWrite   xdmfwrite_
#include <iostream>
//
// C/C++ expect NULL terminated strings. Here is a conversion subroutine.
char * ConvertFortranString( char *FtnName ) {
static char Name[80];
char *np;
memcpy(Name, FtnName, 79 );
Name[79] = '\0';
np = &Name[78];
while( ( np > Name ) && ( *np <= ' ') ) {
       np--;
       }
*np = '\0';
return( Name );
}

//
// C++ will mangle the name based on the argument list. This tells the
// compiler not to mangle the name so we can call it from 'C' (but
// really Fortran in this case)
//
extern "C" {
//

/* 
void
write_xml()
{
    FILE *xmf = 0;
 
    //
    // Open the file and write the XML description of the mesh..
    //
    xmf = fopen("xdmf_animate.xml", "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n", NY+1, NX+1);
    fprintf(xmf, "     <Geometry GeometryType=\"X_Y\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1));
    fprintf(xmf, "        xdmf2d.h5:/X\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1));
    fprintf(xmf, "        xdmf2d.h5:/Y\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    fprintf(xmf, "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
    fprintf(xmf, "        xdmf2d.h5:/Pressure\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "     <Attribute Name=\"VelocityX\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY+1, NX+1);
    fprintf(xmf, "        xdmf2d.h5:/VelocityX\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}
*/

void
XdmfWrite( char *FtnName, int *Iteration,
   int *NumberOfPoints, int *NumberOfTri, XdmfFloat64 *Points,
   XdmfInt32 *Conns, XdmfFloat64 *NodeData,
	   XdmfFloat64 *CellData, int *Outputflags ){
char            *Name;
char            FullName[80];
ostrstream      DataName(FullName, 80);
XdmfDOM         dom;
XdmfRoot        root;
XdmfDomain      domain;
XdmfGrid        grid;
XdmfTime        time;
XdmfTopology    *topology;
XdmfGeometry    *geometry;
XdmfAttribute   nodedata;
XdmfAttribute   celldata;
XdmfArray       *array;
 
 Name = ConvertFortranString( FtnName );
 std::cout << Name;
 // DataName.seekp(0);
 //DataName << Name <<  setfill('0')<< ends;
 //std::cout << FullName;
 root.SetDOM(&dom);
 root.SetVersion(2.0);
 root.Build();
 
 // Domain
 root.Insert(&domain);
 
 // Grid
 grid.SetName("RBC Grid");
 domain.Insert(&grid);
 time.SetTimeType(XDMF_TIME_SINGLE);
 time.SetValue(0.001 * *Iteration);
 grid.Insert(&time);
 
 // Topology
 topology = grid.GetTopology();
 topology->SetTopologyType(XDMF_TRI);
 topology->SetNumberOfElements(*NumberOfTri);
 // Fortran is 1 based while c++ is 0 based so
 // Either subtract 1 from all connections or specify a BaseOffset
 topology->SetBaseOffset(1);
 // If you haven't assigned an XdmfArray, GetConnectivity() will create one.
 array = topology->GetConnectivity();
 array->SetNumberType(XDMF_INT32_TYPE); //XdmfInt32
 array->SetNumberOfElements(*NumberOfTri * 3);
 array->SetValues(0, Conns, *NumberOfTri * 3);
 // C++ string hocus pocus.
 // We're actually building the string in FullName[] but were using streams.
 // the DatasetName will be Demo_00001.h5:/Conns.
 DataName.seekp(0);
 DataName << Name << "_" << setw(9) << setfill('0') << *Iteration << ".h5:/Connectivity" << ends;
 // Where the data will actually be written
 array->SetHeavyDataSetName(FullName);
 
 
 // Geometry
 geometry = grid.GetGeometry();
 geometry->SetGeometryType(XDMF_GEOMETRY_XYZ);
 geometry->SetNumberOfPoints(*NumberOfPoints);
 array = geometry->GetPoints();
 array->SetNumberType(XDMF_FLOAT64_TYPE);
 array->SetValues(0, Points, *NumberOfPoints * 3);
 DataName.seekp(0);
 DataName << Name << "_" << setw(9) << setfill('0') << *Iteration << ".h5:/mesh" << ends;
 array->SetHeavyDataSetName(FullName);
 
 if (Outputflags[2]==1){
   // Node Data
   nodedata.SetName("Node Scalar");
   nodedata.SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_NODE);
   nodedata.SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);
   array = nodedata.GetValues();
   array->SetNumberType(XDMF_FLOAT64_TYPE);
   array->SetNumberOfElements(*NumberOfPoints);
   array->SetValues(0, NodeData, *NumberOfPoints);
   DataName.seekp(0);
   DataName << Name << "_" << setw(9) << setfill('0') << *Iteration << ".h5:/NodeData" << ends;
   array->SetHeavyDataSetName(FullName);
   // Attach and Write
   grid.Insert(&nodedata);
 }
 
 if (Outputflags[3]==1) {
   // Cell Data
   celldata.SetName("Cell Scalar");
   celldata.SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_CELL);
   celldata.SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);
   array = celldata.GetValues();
   array->SetNumberType(XDMF_FLOAT64_TYPE);
   array->SetNumberOfElements(*NumberOfTri);
   array->SetValues(0, CellData, *NumberOfTri);
   DataName.seekp(0);
   DataName << Name << "_" << setw(9) << setfill('0') << *Iteration << ".h5:/CellData" << ends;
   array->SetHeavyDataSetName(FullName);
   // Attach and Write
   grid.Insert(&celldata);
 }
 
 // Attach and Write
 //grid.Insert(&nodedata);
 //grid.Insert(&celldata);
 
// Build is recursive ... it will be called on all of the child nodes.
 // This updates the DOM and writes the HDF5
 root.Build();
 // Write the XML
 DataName.seekp(0);
 DataName << Name << "_" << setw(9) << setfill('0') << *Iteration << ".xmf" << ends;
 dom.Write(FullName);
}
}
