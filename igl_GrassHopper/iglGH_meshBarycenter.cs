﻿using Grasshopper.Kernel;
using System;

namespace igl_Grasshopper
{
    public class IGL_barycenter : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public IGL_barycenter()
          : base("Barycenter", "iBarycenter",
              "compute the barycenter of each triangle of the given mesh.",
              "IGM", "02 | Basic Property")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "input mesh to analysis.", GH_ParamAccess.item);
            //pManager.AddPointParameter("Mesh V", "V", "A list of mesh vertices.", GH_ParamAccess.list);
            //pManager.AddIntegerParameter("Mesh F", "F", "A list of mesh faces.", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Barycenters", "BC", "the barycenters of the input mesh's faces", GH_ParamAccess.list);
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

            // call the cpp function to solve the adjacency list
            var res = IGLRhinoCommon.Utils.getBarycenter(ref mesh);

            // output
            DA.SetDataList(0, res);
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
            get { return new Guid("532f1fab-c27e-4795-99d4-efb55cabc1d4"); }
        }
    }
}
