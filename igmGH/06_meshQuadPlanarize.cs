﻿using Grasshopper.Kernel;
using System;

namespace igmGH
{
    public class IGM_quad_planarize : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public IGM_quad_planarize()
          : base("Planarize Quad Mesh", "igQuadPlanarize",
              "Planarize the quad faces in a quad mesh.",
              "IG-Mesh", "06|Util")
        {
        }

        /// <summary>
        /// icon position in a category
        /// </summary>
        public override GH_Exposure Exposure => GH_Exposure.quarternary;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "Input mesh for analysis.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Iter", "I", "Max iteration for planarization.", GH_ParamAccess.item, 100);
            pManager.AddNumberParameter("Thres", "T", "Threshould to stop the planarization.", GH_ParamAccess.item, 0.005);

            pManager[1].Optional = true;
            pManager[2].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "THe planarized quad mesh.", GH_ParamAccess.item);
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

            int maxIter = 100;
            if (!DA.GetData(1, ref maxIter) || maxIter < 0) { return; }
            double thres = 0.005;
            if (!DA.GetData(2, ref thres) || thres <= 0) { return; }


            // call the cpp function to solve the adjacency list
            IGMRhinoCommon.Utils.planarizeQuadMesh(ref mesh, maxIter, thres);

            // output
            DA.SetData(0, mesh);
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
                return Properties.Resources.meshRandomPtsOnMesh;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("923f1d6e-e4b0-46e2-b913-e78545bfd7ab"); }
        }
    }
}
