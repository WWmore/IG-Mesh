using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace igl_GrassHopper
{
    public class iglGH_perVertFaceNormal : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public iglGH_perVertFaceNormal()
          : base("per_vertex_and_face_normal", "perVertFaceN",
              "compute the mesh normals per vertex and face.",
              "IGL", "Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "input mesh.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddVectorParameter("Vertex Normals", "VN", "vertex-based normals of the mesh", GH_ParamAccess.list);
            pManager.AddVectorParameter("Face Normals", "FN", "face-based normals of the mesh", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Rhino.Geometry.Mesh mesh = new Rhino.Geometry.Mesh();
            if (!DA.GetData(0, ref mesh)) { return; }
            if (!mesh.IsValid) { return; }

            // call the cpp func
            List<Vector3f> VN, FN;
            IGLRhinoCommon.Utils.getPerVertFaceNormal(in mesh, out VN, out FN);

            Grasshopper.DataTree<Vector3f> VNArray = new Grasshopper.DataTree<Vector3f>();
            Grasshopper.DataTree<Vector3f> FNArray = new Grasshopper.DataTree<Vector3f>();

            VNArray.AddRange(VN, new Grasshopper.Kernel.Data.GH_Path(0));
            FNArray.AddRange(FN, new Grasshopper.Kernel.Data.GH_Path(0));

            // set the res on component
            DA.SetDataTree(0, VNArray);
            DA.SetDataTree(1, FNArray);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("EA88C7DD-9DD8-4EAB-89FE-8A79A5D63DF6"); }
        }
    }
}