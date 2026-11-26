
/* Lower56, 26 Nov 2026: exact Gamma72 mixed-shell normal and heavy-profile certificate.
   Run node FILE. No dependencies, solver, external download, or floating point angle tests.
   Mathematical source premises and proof are in lower56_69d_discriminant_torus_20261126.md. */
"use strict";
const SOURCE={"catalogue_gram":[[8,0,0,-4,-4,4,-4,0,-2,4,-4,-2,-4,-2,4,2,-4,4,-4,4,2,-2,-4,-2,-2,-1,0,1,0,-2,2,-1,2,-1,1,2,1,1,0,0,2,1,3,-3,0,1,1,0,2,-1,0,-1,-2,0,0,-1,1,1,-1,1,-1,0,2,1,0,3,1,-1,1,0,-1,-1],[0,8,-4,2,2,2,0,4,-2,0,4,2,4,-4,-2,0,-4,-4,-2,4,-4,4,0,-4,1,-2,2,-2,0,1,-1,-1,1,1,-2,0,-1,0,0,0,-1,2,1,-1,1,0,1,1,1,2,0,-1,1,2,-1,1,0,1,0,1,1,-2,-1,0,-3,0,0,1,-1,2,1,-1],[0,-4,8,-4,2,2,0,-4,4,2,-2,2,-4,2,4,-2,2,2,4,-2,4,-2,0,2,0,0,-2,3,-2,-1,-1,0,-2,-1,1,-3,1,0,0,-1,1,-1,-1,0,-1,1,-2,0,0,-2,2,1,-1,0,-1,-2,0,0,0,-2,-1,1,2,-2,2,0,1,-1,1,0,-2,1],[-4,2,-4,8,0,-2,0,0,-2,-4,4,-2,4,-2,-4,-2,0,-4,0,-2,-4,4,0,-2,1,1,-1,-2,1,0,1,2,-1,0,-1,0,0,0,-2,2,-1,0,-2,2,-1,-2,1,1,-1,2,-3,2,1,-1,1,2,-2,-2,1,-1,2,-1,-4,1,-1,-2,-2,1,-3,0,1,0],[-4,2,2,0,8,0,4,0,4,-2,2,4,2,2,-2,0,0,-4,4,0,-2,0,2,2,2,-1,1,-1,-2,2,-2,0,-1,3,-2,-2,-1,0,1,-2,-1,1,-2,1,1,0,-2,0,0,0,2,-1,2,2,0,0,1,2,-1,0,0,1,0,-2,-1,-1,0,1,0,0,-1,1],[4,2,2,-2,0,8,-4,0,0,4,-2,0,-2,-4,2,-2,-4,2,-2,4,2,2,-4,-4,0,-2,0,1,-2,-2,0,-1,0,0,0,-1,0,1,0,0,1,1,1,-2,-1,0,0,0,2,-1,1,0,-2,2,-2,-1,0,2,-1,-1,-1,-1,1,-1,-1,2,0,0,0,1,-2,-2],[-4,0,0,0,4,-4,8,2,4,-2,0,4,2,4,-2,2,4,-2,2,-2,-2,-2,2,4,0,1,1,-1,0,2,-2,1,-1,2,0,-1,0,-1,1,-2,-1,-1,-1,2,1,1,0,0,-2,1,1,-1,2,0,2,2,1,1,0,1,1,1,0,-1,1,-2,0,1,0,0,1,2],[0,4,-4,0,0,0,2,8,0,2,0,0,4,0,-2,0,0,-2,-4,4,-2,0,2,0,1,-1,2,-2,0,1,-2,-2,0,1,-1,0,-2,-1,2,1,-2,1,1,-1,1,0,1,0,1,1,0,-2,0,1,-1,2,0,2,-1,0,0,-1,1,1,-2,0,-1,1,0,0,2,0],[-2,-2,4,-2,4,0,4,0,8,2,-2,4,0,4,2,-2,2,0,4,-2,2,-2,2,2,-1,0,0,2,-1,0,-1,0,-2,0,1,-2,1,-1,0,-2,0,-1,-1,0,-1,1,-1,1,-2,-1,2,1,1,0,1,0,2,1,0,0,1,1,1,-3,1,-1,1,-1,0,0,0,2],[4,0,2,-4,-2,4,-2,2,2,8,-4,0,-2,0,4,-2,-2,2,-2,2,2,-2,-2,-2,-1,-1,0,2,-2,-2,-1,-2,-1,-2,2,-1,1,-1,1,0,1,0,1,-3,0,2,1,0,1,-1,1,0,-3,0,-2,-1,0,2,0,-1,0,-1,3,-1,0,1,0,-2,1,1,0,-1],[-4,4,-2,4,2,-2,0,0,-2,-4,8,2,4,-2,-2,-2,0,-4,2,0,-2,4,2,-2,1,0,0,-1,1,1,0,1,0,0,-2,0,-1,0,-2,1,-1,0,-1,1,0,-1,1,1,-1,2,-1,1,2,0,0,1,-1,-2,2,1,1,-1,-3,0,-1,-2,0,1,-1,1,2,0],[-2,2,2,-2,4,0,4,0,4,0,2,8,0,0,0,0,0,0,4,0,0,0,2,0,-1,-1,2,1,0,1,-1,0,0,1,-1,-2,0,0,0,-2,1,-1,0,-1,0,1,0,1,-2,0,3,0,2,1,1,0,2,1,0,2,0,0,0,-2,1,-1,2,-1,0,1,1,1],[-4,4,-4,4,2,-2,2,4,0,-2,4,0,8,0,-4,-2,0,-4,0,0,-4,4,2,0,1,-1,1,-2,0,1,-1,0,-1,0,-1,0,-2,-1,0,1,-1,1,-1,1,1,-1,1,1,-1,1,-1,0,1,0,0,2,-1,-1,1,0,2,-1,-2,0,-1,-1,-1,1,-1,1,2,1],[-2,-4,2,-2,2,-4,4,0,4,0,-2,0,0,8,0,0,4,0,2,-2,2,-4,2,4,0,2,-1,1,-1,1,-1,1,-1,1,1,0,1,-2,1,-2,-1,-1,-1,1,0,1,-1,-1,-1,0,0,0,0,-1,1,1,1,1,0,0,1,2,1,-2,1,-1,0,0,1,-1,0,1],[4,-2,4,-4,-2,2,-2,-2,2,4,-2,0,-4,0,8,0,0,2,0,0,4,-2,-2,0,-2,1,-2,4,0,-1,0,-1,-1,-3,3,0,2,-1,-2,-1,0,-1,1,-1,-1,3,0,-1,0,0,0,2,-1,0,-1,-2,0,-1,2,0,0,-1,2,-1,0,0,1,-1,1,2,-1,-1],[2,0,-2,-2,0,-2,2,0,-2,-2,-2,0,-2,0,0,8,0,0,-2,0,-2,-2,0,2,-1,0,2,-1,2,1,1,-1,3,1,0,2,0,2,1,-2,0,0,2,0,2,1,0,-1,0,0,1,-2,2,0,2,-1,2,0,-1,2,-1,2,1,2,0,0,1,0,1,0,0,0],[-4,-4,2,0,0,-4,4,0,2,-2,0,0,0,4,0,0,8,0,2,-4,2,-2,2,4,0,3,-2,1,1,1,-1,2,-1,0,1,-1,1,-1,0,0,-2,-3,-2,3,-1,1,0,-1,-2,1,-1,1,1,-1,1,2,0,-1,1,-1,1,1,0,0,2,-3,-1,1,0,0,1,1],[4,-4,2,-4,-4,2,-2,-2,0,2,-4,0,-4,0,2,0,0,8,0,0,4,-2,-2,0,-3,0,0,2,1,-2,2,0,1,-1,2,1,1,1,0,0,3,-2,2,-2,-1,0,0,0,-1,-2,1,0,-1,-1,1,-1,1,0,0,1,-1,1,1,0,3,2,2,-2,1,-1,-1,0],[-4,-2,4,0,4,-2,2,-4,4,-2,2,4,0,2,0,-2,2,0,8,-4,0,0,2,2,-1,0,-1,2,0,0,0,1,-1,0,0,-2,1,0,-1,-1,1,-2,-2,0,-1,0,-1,1,-3,-1,1,2,2,-1,1,-1,1,-1,1,0,1,1,-1,-2,2,-2,2,-2,-1,0,0,2],[4,4,-2,-2,0,4,-2,4,-2,2,0,0,0,-2,0,0,-4,0,-4,8,0,0,-2,-2,1,-1,1,-1,-1,0,-1,-1,1,2,-1,1,-1,0,1,0,-1,2,2,-2,0,0,1,-1,3,1,0,-2,-1,2,-2,1,0,3,-1,1,-1,-1,1,0,-3,2,0,2,0,0,0,-2],[2,-4,4,-4,-2,2,-2,-2,2,2,-2,0,-4,2,4,-2,2,4,0,0,8,-2,0,0,-1,1,-1,3,0,0,0,0,0,-1,1,0,1,-1,-1,-1,0,-1,1,0,-2,1,-1,-1,0,-1,1,1,-1,1,-1,-1,1,0,0,0,-1,0,1,-2,1,1,1,0,2,0,-1,-1],[-2,4,-2,4,0,2,-2,0,-2,-2,4,0,4,-4,-2,-2,-2,-2,0,0,-2,8,0,-4,0,-2,0,0,0,-1,0,0,0,-1,-1,-1,-1,1,-2,0,0,1,0,0,0,-2,1,1,-1,0,-1,2,0,0,-1,0,-1,-2,1,-1,1,-1,-3,-1,-1,0,0,0,-1,2,1,-1],[-4,0,0,0,2,-4,2,2,2,-2,2,2,2,2,-2,0,2,-2,2,-2,0,0,8,2,1,-1,2,-1,1,2,-1,-2,0,0,-2,-1,-2,0,1,0,-1,1,0,0,1,-1,-2,1,-1,-1,2,-1,2,0,0,-1,1,-1,-1,0,-1,1,0,0,0,0,1,-1,1,-1,2,2],[-2,-4,2,-2,2,-4,4,0,2,-2,-2,0,0,4,0,2,4,0,2,-2,0,-4,2,8,1,1,-1,0,-1,2,-2,0,-2,1,0,-1,-1,-1,1,0,-1,0,-2,2,1,1,-2,-2,0,-1,0,-1,0,0,0,0,-1,0,-1,-1,-1,1,1,1,1,0,-1,1,1,-1,-1,2],[-2,1,0,1,2,0,0,1,-1,-1,1,-1,1,0,-2,-1,0,-3,-1,1,-1,0,1,1,8,0,0,-4,-4,4,-4,0,-2,4,-4,-2,-4,-2,4,2,-4,4,-4,4,2,-2,-4,-2,2,1,0,-1,0,2,-2,1,-2,1,-1,-2,-1,-1,0,0,-2,-1,-3,3,0,-1,-1,0],[-1,-2,0,1,-1,-2,1,-1,0,-1,0,-1,-1,2,1,0,3,0,0,-1,1,-2,-1,1,0,8,-4,2,2,2,0,4,-2,0,4,2,4,-4,-2,0,-4,-4,-2,4,-4,4,0,-4,-1,2,-2,2,0,-1,1,1,-1,-1,2,0,1,0,0,0,1,-2,-1,1,-1,0,-1,-1],[0,2,-2,-1,1,0,1,2,0,0,0,2,1,-1,-2,2,-2,0,-1,1,-1,0,2,-1,0,-4,8,-4,2,2,0,-4,4,2,-2,2,-4,2,4,-2,2,2,4,-2,4,-2,0,2,0,0,2,-3,2,1,1,0,2,1,-1,3,-1,0,0,1,-1,1,1,0,1,-1,2,0],[1,-2,3,-2,-1,1,-1,-2,2,2,-1,1,-2,1,4,-1,1,2,2,-1,3,0,-1,0,-4,2,-4,8,0,-2,0,0,-2,-4,4,-2,4,-2,-4,-2,0,-4,0,-2,-4,4,0,-2,-1,-1,1,2,-1,0,-1,-2,1,0,1,0,0,0,2,-2,1,0,2,-2,1,2,-1,-1],[0,0,-2,1,-2,-2,0,0,-1,-2,1,0,0,-1,0,2,1,1,0,-1,0,0,1,-1,-4,2,2,0,8,0,4,0,4,-2,2,4,2,2,-2,0,0,-4,4,0,-2,0,2,2,-2,1,-1,1,2,-2,2,0,1,-3,2,2,1,0,-1,2,1,-1,2,-1,-1,0,2,0],[-2,1,-1,0,2,-2,2,1,0,-2,1,1,1,1,-1,1,1,-2,0,0,0,-1,2,2,4,2,2,-2,0,8,-4,0,0,4,-2,0,-2,-4,2,-2,-4,2,-2,4,2,2,-4,-4,0,2,0,-1,2,2,0,1,0,0,0,1,0,-1,0,0,-1,-1,-1,2,1,0,0,0],[2,-1,-1,1,-2,0,-2,-2,-1,-1,0,-1,-1,-1,0,1,-1,2,0,-1,0,0,-1,-2,-4,0,0,0,4,-4,8,2,4,-2,0,4,2,4,-2,2,4,-2,2,-2,-2,-2,2,4,0,-1,-1,1,0,-2,2,-1,1,-2,0,1,0,1,-1,2,1,1,1,-2,-1,-1,0,0],[-1,-1,0,2,0,-1,1,-2,0,-2,1,0,0,1,-1,-1,2,0,1,-1,0,0,-2,0,0,4,-4,0,0,0,2,8,0,2,0,0,4,0,-2,0,0,-2,-4,4,-2,0,2,0,-1,1,-2,2,0,-1,2,2,0,-1,1,0,2,1,-2,-1,2,-1,-1,1,-1,0,-1,0],[2,1,-2,-1,-1,0,-1,0,-2,-1,0,0,-1,-1,-1,3,-1,1,-1,1,0,0,0,-2,-2,-2,4,-2,4,0,4,0,8,2,-2,4,0,4,2,-2,2,0,4,-2,2,-2,2,2,1,0,0,-2,1,0,1,0,2,0,-1,2,-1,1,0,2,0,1,1,0,1,-1,1,-1],[-1,1,-1,0,3,0,2,1,0,-2,0,1,0,1,-3,1,0,-1,0,2,-1,-1,0,1,4,0,2,-4,-2,4,-2,2,2,8,-4,0,-2,0,4,-2,-2,2,-2,2,2,-2,-2,-2,1,1,0,-2,2,2,1,2,1,2,-2,1,-1,1,-1,0,-1,0,-1,3,0,-2,-1,0],[1,-2,1,-1,-2,0,0,-1,1,2,-2,-1,-1,1,3,0,1,2,0,-1,1,-1,-2,0,-4,4,-2,4,2,-2,0,0,-2,-4,8,2,4,-2,-2,-2,0,-4,2,0,-2,4,2,-2,-1,0,0,1,-1,-1,0,-1,0,0,2,0,1,0,2,-1,1,0,1,-1,0,1,-1,-1],[2,0,-3,0,-2,-1,-1,0,-2,-1,0,-2,0,0,0,2,-1,1,-2,1,0,-1,-1,-1,-2,2,2,-2,4,0,4,0,4,0,2,8,0,0,0,0,0,0,4,0,0,0,2,0,1,1,-2,-1,0,-1,1,0,0,-1,1,2,0,0,0,2,-1,1,0,1,0,-1,0,-1],[1,-1,1,0,-1,0,0,-2,1,1,-1,0,-2,1,2,0,1,1,1,-1,1,-1,-2,-1,-4,4,-4,4,2,-2,2,4,0,-2,4,0,8,0,-4,-2,0,-4,0,0,-4,4,2,0,-1,1,-1,2,0,-1,1,0,1,0,1,0,2,1,0,-1,1,-1,1,-1,-1,1,-1,-1],[1,0,0,0,0,1,-1,-1,-1,-1,0,0,-1,-2,-1,2,-1,1,0,0,-1,1,0,-1,-2,-4,2,-2,2,-4,4,0,4,0,-2,0,0,8,0,0,4,0,2,-2,2,-4,2,4,0,-2,1,-1,1,-1,1,-1,1,-1,-1,0,-1,2,-1,2,1,1,1,-1,0,-1,1,1],[0,0,0,-2,1,0,1,2,0,1,-2,0,0,1,-2,1,0,0,-1,1,-1,-2,1,1,4,-2,4,-4,-2,2,-2,-2,2,4,-2,0,-4,0,8,0,0,2,0,0,4,-2,-2,0,2,-1,2,-4,0,1,0,1,1,3,-3,0,-2,1,2,1,0,1,-1,1,1,-3,0,1],[0,0,-1,2,-2,0,-2,1,-2,0,1,-2,1,-2,-1,-2,0,0,-1,0,-1,0,0,0,2,0,-2,-2,0,-2,2,0,-2,-2,-2,0,-2,0,0,8,0,0,-2,0,-2,-2,0,2,1,0,-2,1,-2,-1,-1,1,-3,-1,0,-2,0,-2,-1,2,0,0,-2,0,-2,-1,0,1],[2,-1,1,-1,-1,1,-1,-2,0,1,-1,1,-1,-1,0,0,-2,3,1,-1,0,0,-1,-1,-4,-4,2,0,0,-4,4,0,2,-2,0,0,0,4,0,0,8,0,2,-4,2,-2,2,4,0,-3,2,-1,-1,-1,1,-2,1,0,-1,1,-1,1,0,0,2,3,2,-3,1,-1,0,1],[1,2,-1,0,1,1,-1,1,-1,0,0,-1,1,-1,-1,0,-3,-2,-2,2,-1,1,1,0,4,-4,2,-4,-4,2,-2,-2,0,2,-4,0,-4,0,2,0,0,8,0,0,4,-2,-2,0,3,0,0,-2,-1,2,-2,0,-1,1,-2,-1,-1,-1,0,0,-3,2,-2,2,1,0,0,0],[3,1,-1,-2,-2,1,-1,1,-1,1,-1,0,-1,-1,1,2,-2,2,-2,2,1,0,0,-2,-4,-2,4,0,4,-2,2,-4,4,-2,2,4,0,2,0,-2,2,0,8,-4,0,0,2,2,1,0,1,-2,0,0,0,-1,1,0,0,2,-1,0,1,1,-1,2,2,0,1,0,1,-1],[-3,-1,0,2,1,-2,2,-1,0,-3,1,-1,1,1,-1,0,3,-2,0,-2,0,0,0,2,4,4,-2,-2,0,4,-2,4,-2,2,0,0,0,-2,0,0,-4,0,-4,8,0,0,-2,-2,-1,1,-1,1,1,0,1,1,-1,-2,1,-1,1,0,-1,0,1,-2,-2,2,0,0,-1,1],[0,1,-1,-1,1,-1,1,1,-1,0,0,0,1,0,-1,2,-1,-1,-1,0,-2,0,1,1,2,-4,4,-4,-2,2,-2,-2,2,2,-2,0,-4,2,4,-2,2,4,0,0,8,-2,0,0,1,-1,1,-3,0,0,0,0,0,1,-1,0,-1,1,1,1,0,1,-1,0,2,-1,1,1],[1,0,1,-2,0,0,1,0,1,2,-1,1,-1,1,3,1,1,0,0,0,1,-2,-1,1,-2,4,-2,4,0,2,-2,0,-2,-2,4,0,4,-4,-2,-2,-2,-2,0,0,-2,8,0,-4,0,2,0,0,0,1,0,0,0,1,1,1,1,-1,2,0,0,-1,0,0,0,2,-1,-1],[1,1,-2,1,-2,0,0,1,-1,1,1,0,1,-1,0,0,0,0,-1,1,-1,1,-2,-2,-4,0,0,0,2,-4,2,2,2,-2,2,2,2,2,-2,0,2,-2,2,-2,0,0,8,2,-1,1,-2,1,-1,-2,1,2,0,0,2,1,2,0,-1,0,1,-1,0,0,-1,1,2,-1],[0,1,0,1,0,0,0,0,1,0,1,1,1,-1,-1,-1,-1,0,1,-1,-1,1,1,-2,-2,-4,2,-2,2,-4,4,0,2,-2,-2,0,0,4,0,2,4,0,2,-2,0,-4,2,8,-1,-1,1,0,1,-2,2,0,2,-1,0,1,1,1,-1,0,1,0,2,-2,-1,-1,2,2],[2,1,0,-1,0,2,-2,1,-2,1,-1,-2,-1,-1,0,0,-2,-1,-3,3,0,-1,-1,0,2,-1,0,-1,-2,0,0,-1,1,1,-1,1,-1,0,2,1,0,3,1,-1,1,0,-1,-1,8,0,0,-4,-4,4,-4,0,-2,4,-4,-2,-4,-2,4,2,-4,4,-4,4,2,-2,-4,-2],[-1,2,-2,2,0,-1,1,1,-1,-1,2,0,1,0,0,0,1,-2,-1,1,-1,0,-1,-1,1,2,0,-1,1,2,-1,1,0,1,0,1,1,-2,-1,0,-3,0,0,1,-1,2,1,-1,0,8,-4,2,2,2,0,4,-2,0,4,2,4,-4,-2,0,-4,-4,-2,4,-4,4,0,-4],[0,0,2,-3,2,1,1,0,2,1,-1,3,-1,0,0,1,-1,1,1,0,1,-1,2,0,0,-2,2,1,-1,0,-1,-2,0,0,0,-2,-1,1,2,-2,2,0,1,-1,1,0,-2,1,0,-4,8,-4,2,2,0,-4,4,2,-2,2,-4,2,4,-2,2,2,4,-2,4,-2,0,2],[-1,-1,1,2,-1,0,-1,-2,1,0,1,0,0,0,2,-2,1,0,2,-2,1,2,-1,-1,-1,2,-3,2,1,-1,1,2,-2,-2,1,-1,2,-1,-4,1,-1,-2,-2,1,-3,0,1,0,-4,2,-4,8,0,-2,0,0,-2,-4,4,-2,4,-2,-4,-2,0,-4,0,-2,-4,4,0,-2],[-2,1,-1,1,2,-2,2,0,1,-3,2,2,1,0,-1,2,1,-1,2,-1,-1,0,2,0,0,0,2,-1,2,2,0,0,1,2,-1,0,0,1,0,-2,-1,-1,0,1,0,0,-1,1,-4,2,2,0,8,0,4,0,4,-2,2,4,2,2,-2,0,0,-4,4,0,-2,0,2,2],[0,2,0,-1,2,2,0,1,0,0,0,1,0,-1,0,0,-1,-1,-1,2,1,0,0,0,2,-1,1,0,-2,2,-2,-1,0,2,-1,-1,-1,-1,1,-1,-1,2,0,0,0,1,-2,-2,4,2,2,-2,0,8,-4,0,0,4,-2,0,-2,-4,2,-2,-4,2,-2,4,2,2,-4,-4],[0,-1,-1,1,0,-2,2,-1,1,-2,0,1,0,1,-1,2,1,1,1,-2,-1,-1,0,0,-2,1,1,-1,2,0,2,2,1,1,0,1,1,1,0,-1,1,-2,0,1,0,0,1,2,-4,0,0,0,4,-4,8,2,4,-2,0,4,2,4,-2,2,4,-2,2,-2,-2,-2,2,4],[-1,1,-2,2,0,-1,2,2,0,-1,1,0,2,1,-2,-1,2,-1,-1,1,-1,0,-1,0,1,1,0,-2,0,1,-1,2,0,2,-1,0,0,-1,1,1,-2,0,-1,1,0,0,2,0,0,4,-4,0,0,0,2,8,0,2,0,0,4,0,-2,0,0,-2,-4,4,-2,0,2,0],[1,0,0,-2,1,0,1,0,2,0,-1,2,-1,1,0,2,0,1,1,0,1,-1,1,-1,-2,-1,2,1,1,0,1,0,2,1,0,0,1,1,1,-3,1,-1,1,-1,0,0,0,2,-2,-2,4,-2,4,0,4,0,8,2,-2,4,0,4,2,-2,2,0,4,-2,2,-2,2,2],[1,1,0,-2,2,2,1,2,1,2,-2,1,-1,1,-1,0,-1,0,-1,3,0,-2,-1,0,1,-1,1,0,-3,0,-2,-1,0,2,0,-1,0,-1,3,-1,0,1,0,-2,1,1,0,-1,4,0,2,-4,-2,4,-2,2,2,8,-4,0,-2,0,4,-2,-2,2,-2,2,2,-2,-2,-2],[-1,0,0,1,-1,-1,0,-1,0,0,2,0,1,0,2,-1,1,0,1,-1,0,1,-1,-1,-1,2,-1,1,2,0,0,1,-1,-2,2,1,1,-1,-3,0,-1,-2,0,1,-1,1,2,0,-4,4,-2,4,2,-2,0,0,-2,-4,8,2,4,-2,-2,-2,0,-4,2,0,-2,4,2,-2],[1,1,-2,-1,0,-1,1,0,0,-1,1,2,0,0,0,2,-1,1,0,1,0,-1,0,-1,-2,0,3,0,2,1,1,0,2,1,0,2,0,0,0,-2,1,-1,2,-1,0,1,1,1,-2,2,2,-2,4,0,4,0,4,0,2,8,0,0,0,0,0,0,4,0,0,0,2,0],[-1,1,-1,2,0,-1,1,0,1,0,1,0,2,1,0,-1,1,-1,1,-1,-1,1,-1,-1,-1,1,-1,0,1,0,0,2,-1,-1,1,0,2,-1,-2,0,-1,-1,-1,1,-1,1,2,1,-4,4,-4,4,2,-2,2,4,0,-2,4,0,8,0,-4,-2,0,-4,0,0,-4,4,2,0],[0,-2,1,-1,1,-1,1,-1,1,-1,-1,0,-1,2,-1,2,1,1,1,-1,0,-1,1,1,-1,0,0,0,0,-1,1,1,1,1,0,0,1,2,1,-2,1,-1,0,0,1,-1,0,1,-2,-4,2,-2,2,-4,4,0,4,0,-2,0,0,8,0,0,4,0,2,-2,2,-4,2,4],[2,-1,2,-4,0,1,0,1,1,3,-3,0,-2,1,2,1,0,1,-1,1,1,-3,0,1,0,0,0,2,-1,0,-1,-2,0,-1,2,0,0,-1,2,-1,0,0,1,-1,1,2,-1,-1,4,-2,4,-4,-2,2,-2,-2,2,4,-2,0,-4,0,8,0,0,2,0,0,4,-2,-2,0],[1,0,-2,1,-2,-1,-1,1,-3,-1,0,-2,0,-2,-1,2,0,0,-2,0,-2,-1,0,1,0,0,1,-2,2,0,2,-1,2,0,-1,2,-1,2,1,2,0,0,1,0,1,0,0,0,2,0,-2,-2,0,-2,2,0,-2,-2,-2,0,-2,0,0,8,0,0,-2,0,-2,-2,0,2],[0,-3,2,-1,-1,-1,1,-2,1,0,-1,1,-1,1,0,0,2,3,2,-3,1,-1,0,1,-2,1,-1,1,1,-1,1,2,0,-1,1,-1,1,1,0,0,2,-3,-1,1,0,0,1,1,-4,-4,2,0,0,-4,4,0,2,-2,0,0,0,4,0,0,8,0,2,-4,2,-2,2,4],[3,0,0,-2,-1,2,-2,0,-1,1,-2,-1,-1,-1,0,0,-3,2,-2,2,1,0,0,0,-1,-2,1,0,-1,-1,1,-1,1,0,0,1,-1,1,1,0,3,2,2,-2,1,-1,-1,0,4,-4,2,-4,-4,2,-2,-2,0,2,-4,0,-4,0,2,0,0,8,0,0,4,-2,-2,0],[1,0,1,-2,0,0,0,-1,1,0,0,2,-1,0,1,1,-1,2,2,0,1,0,1,-1,-3,-1,1,2,2,-1,1,-1,1,-1,1,0,1,1,-1,-2,2,-2,2,-2,-1,0,0,2,-4,-2,4,0,4,-2,2,-4,4,-2,2,4,0,2,0,-2,2,0,8,-4,0,0,2,2],[-1,1,-1,1,1,0,1,1,-1,-2,1,-1,1,0,-1,0,1,-2,-2,2,0,0,-1,1,3,1,0,-2,-1,2,-2,1,0,3,-1,1,-1,-1,1,0,-3,2,0,2,0,0,0,-2,4,4,-2,-2,0,4,-2,4,-2,2,0,0,0,-2,0,0,-4,0,-4,8,0,0,-2,-2],[1,-1,1,-3,0,0,0,0,0,1,-1,0,-1,1,1,1,0,1,-1,0,2,-1,1,1,0,-1,1,1,-1,1,-1,-1,1,0,0,0,-1,0,1,-2,1,1,1,0,2,0,-1,-1,2,-4,4,-4,-2,2,-2,-2,2,2,-2,0,-4,2,4,-2,2,4,0,0,8,-2,0,0],[0,2,0,0,0,1,0,0,0,1,1,1,1,-1,2,0,0,-1,0,0,0,2,-1,-1,-1,0,-1,2,0,0,-1,0,-1,-2,1,-1,1,-1,-3,-1,-1,0,0,0,-1,2,1,-1,-2,4,-2,4,0,2,-2,0,-2,-2,4,0,4,-4,-2,-2,-2,-2,0,0,-2,8,0,-4],[-1,1,-2,1,-1,-2,1,2,0,0,2,1,2,0,-1,0,1,-1,0,0,-1,1,2,-1,-1,-1,2,-1,2,0,0,-1,1,-1,-1,0,-1,1,0,0,0,0,1,-1,1,-1,2,2,-4,0,0,0,2,-4,2,2,2,-2,2,2,2,2,-2,0,2,-2,2,-2,0,0,8,2],[-1,-1,1,0,1,-2,2,0,2,-1,0,1,1,1,-1,0,1,0,2,-2,-1,-1,2,2,0,-1,0,-1,0,0,0,0,-1,0,-1,-1,-1,1,1,1,1,0,-1,1,1,-1,-1,2,-2,-4,2,-2,2,-4,4,0,2,-2,-2,0,0,4,0,2,4,0,2,-2,0,-4,2,8]],"normals":[[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],"source_blobs":{"gamma72_gram.npy":"5857c4e610c14b6e6e56e734e5102546f356c1dc","gamma72_d4_k4.npy":"51e91c8f489ef600b4e6b8e054654c2ea06650df"},"sha256":{"gamma72_gram.npy":"ddd7bd1d4e4b2005181efa54ad5bf95fbde287a4dfa7261f412ffb3c2b384535","gamma72_d4_k4.npy":"e29723a85fb84bcb37165ef75aaa610b9566ecb255dc54c25c1adc8b37159789"}};
const makeA3=function(){
const mul=(a,b)=>a.map((_,i)=>a.slice(0,i+1).reduce((z,v,j)=>z+v*b[i-j],0n));
const pow=(a,e)=>{let r=a.map((_,i)=>i?0n:1n);for(let i=0;i<e;i++)r=mul(a,r);return r};
const K=15;
const one=Array.from({length:K},(_,i)=>i?0n:1n);
const E4=one.map((_,i)=>i?240n*Array.from({length:i},(_,j)=>j+1).filter(j=>i%j===0).reduce((s,j)=>s+BigInt(j**3),0n):1n);
let Del=one.map((_,i)=>i===1?1n:0n);for(let j=1;j<K;j++){const a=one.map((x,i)=>i===j?-1n:x);Del=mul(Del,pow(a,24));}
const bases=Array.from({length:4},(_,i)=>mul(pow(E4,9-3*i),pow(Del,i)));
let theta=bases[0];let tcoef=[1n];
for(let i=1;i<4;i++){let z=-theta[i];tcoef.push(z);theta=theta.map((x,j)=>x+z*bases[i][j]);}
function q(h){return h.reduce((s,x)=>s+x*x,0)}
function minrep(h){let m=Infinity,vr=[];for(let a=-4;a<=4;a++)for(let b=-4;b<=4;b++)for(let c=-4;c<=4;c++)if((a+b+c)%2===0){let t=h.map((v,i)=>v+8*[a,b,c][i]);let z=q(t);if(z<m){m=z;vr=[t]}else if(z===m)vr.push(t);}
return [m,vr];}
const red=new Map();function cached(h){let k=h.join();if(!red.has(k))red.set(k,minrep(h)); return red.get(k);}
function sort(h){return h.map(Math.abs).sort((a,b)=>b-a);}
const shell=[8,10,12];
let P=[], M=new Map(), U=[];
for(let s of shell)for(let a=-s;a<=s;a++)for(let b=-s;b<=s;b++)for(let c=-s;c<=s;c++){
let h=[a,b,c];if(((a^b)|(a^c))&1)continue;
if(Math.abs(a)+Math.abs(b)>s||Math.abs(a)+Math.abs(c)>s||Math.abs(b)+Math.abs(c)>s)continue;
let r=16*s-q(h); if(r<=0)continue;
let [m,reps]=cached(h);if(r+m<128)continue;
let can=reps.map(sort).sort((a,b)=>a[0]-b[0]||a[1]-b[1]||a[2]-b[2])[0];
let key=r+":"+can.join();
if(!M.has(key)){M.set(key,U.length);U.push({r,can,s0:(r+m)/16,nr:0});}
let i=M.get(key);P.push({s,h,r,i,m,origin:q(h)===m});if(q(h)===m)U[i].nr++;
}
let C=P.filter(z=>z.origin&&z.r<=128&&z.h.join()===cached(z.h)[1].map(x=>x.join()).sort()[0]);
let edge=C.map(a=>C.map(b=>{
let [m]=cached(a.h.map((x,i)=>x-b.h[i]));const x=a.r+b.r-128+m;
return +(x<=0||x*x<=a.r*b.r);
}));
const ps=[];for(let i=0;i<=5;i++)for(let j=0;j<=i;j++)for(let k=0;k<=j;k++){if(i+j+k<=5)ps.push([i,j,k]);}
const roots=[];for(let i=0;i<3;i++)for(let j=i+1;j<3;j++)for(let a of [-8,8])for(let b of [-8,8]){let r=[0,0,0];r[i]=a;r[j]=b;roots.push(r);}
const A=[],B=[],raw=[];
function f(h,p){return h.reduce((x,v,j)=>x*v**(2*p[j]),1);}
for(let s of shell)for(let p of ps){
let row=U.map(()=>0);for(let z of P)if(z.s===s)row[z.i]+=f(z.h,p);
let l=p.reduce((a,b)=>a+b,0);let num=theta[s/2]*(16n*BigInt(s))**BigInt(l),den=1n;
for(let a of p)for(let j=1;j<=a;j++)num*=BigInt(2*j-1);
for(let j=0;j<l;j++)den*=BigInt(72+2*j);
if(s===8)for(let h of roots)num-=den*BigInt(f(h,p));
A.push(row);B.push(Number(num)/Number(den));raw.push({s,p,num:String(num),den:String(den)});
}
return {U,P,C,edge,A,B,raw,roots,theta:theta.map(String),tcoef:tcoef.map(String)};
};
const makeD4=function(thetaAll){
const n=4,S=[8,10,12];const theta=thetaAll;
const q=h=>h.reduce((s,x)=>s+x*x,0);
const H=[[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]];
const hfun=h=>H.map(r=>r.reduce((t,x,i)=>t+x*h[i],0)/2);
const sort=h=>h.map(Math.abs).sort((a,b)=>b-a);
const lex=(a,b)=>{for(let i=0;i<n;i++)if(a[i]!==b[i])return a[i]-b[i];return 0};
function minrep(h){
 let v=h.map(x=>Math.floor((x+4)/8));let r=h.map((x,i)=>x-8*v[i]);let k=v.reduce((a,b)=>a+b,0)&1;
 if(k){let max=-1;for(let i=0;i<n;i++)if(max<0||Math.abs(r[i])>Math.abs(r[max]))max=i;
 r[max]+=r[max]>0?-8:8;}
 return r;
}
function minQ(h){return q(minrep(h));}
function key(h){let v=h.map(x=>Math.floor(x/8));let r=h.map((x,i)=>x-8*v[i]);let k=v.reduce((a,b)=>a+b,0)&1;return r.join()+":"+k;}
let canonical=new Map();function orbit(h){
 let b=sort(h);let id=b.join();if(canonical.has(id))return canonical.get(id);
 let c=sort(hfun(b)); let d=sort(hfun([-b[0],...b.slice(1)]));let r=[b,c,d].sort(lex)[0];
 canonical.set(id,r);return r;
}
let low=new Map(),P=[],M=new Map(),U=[];
function elk(s,h,r,m){
 if(r+m<128)return;
 let rep=minrep(h),ca=orbit(rep);
 let id=r+":"+ca.join();if(!M.has(id)){M.set(id,U.length);U.push({r,can:ca,s0:(r+m)/16});}
 let i=M.get(id);P.push({s,h,r,i,m,origin:q(h)===m});
 if(q(h)===m){let k=r+":"+key(h),prev=low.get(k);if(!prev||lex(h,prev.h)<0)low.set(k,{s,h,r,i,m,origin:true});}
}
let R0=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let a of [-8,8])for(let b of [-8,8]){let r=Array(n).fill(0);r[i]=a;r[j]=b;R0.push([8,r]);}
const p0=Array(n).fill(0);
for(let s of S){
function rec(k,a,b,parity){
 if(k===n){
  if(a>16*s)return;
  if(a===16*s){if(p0.every(x=>x%8===0)&&p0.reduce((x,y)=>x+y,0)%16===0)R0.push([s,p0.slice()]);return;}
  let h=p0.slice(),r=16*s-a,m=minQ(h);elk(s,h,r,m);return;
 }
 for(let t=-s;t<=s;t++)if((k===0||((t&1)===parity))&&Math.abs(t)+b<=s){
  p0[k]=t;rec(k+1,a+t*t,Math.max(b,Math.abs(t)),t&1);
 }
}
rec(0,0,0,0);
}
let C=[...low.values()].filter(z=>z.r<=128);
let ps=[];
for(let t=0;t<=5;t++)for(let i=0;i<=t;i++)for(let j=0;j<=i;j++)for(let k=0;k<=j;k++){let ell=t-i-j-k;if(ell>=0&&ell<=k)ps.push([i,j,k,ell]);}
const roots=R0.map(x=>x[1]);const A=[],B=[],raw=[];
const f=(h,p)=>h.reduce((x,v,j)=>x*v**(2*p[j]),1);
for(let s of S)for(let p of ps){
 let row=U.map(()=>0);for(let z of P)if(z.s===s)row[z.i]+=f(z.h,p);
 let l=p.reduce((a,b)=>a+b,0);let num=BigInt(theta[s/2])*(16n*BigInt(s))**BigInt(l),den=1n;
 for(let a of p)for(let j=1;j<=a;j++)num*=BigInt(2*j-1);
 for(let j=0;j<l;j++)den*=BigInt(72+2*j);
 if(s===8)for(let h of roots)num-=den*BigInt(f(h,p));
 A.push(row);B.push(Number(num)/Number(den));raw.push({s,p,num:String(num),den:String(den)});
}
return {U,P,C,A,B,raw,roots,theta,low:R0,minQ, D4:true};
};
const harmA3=function(model){
const factor=[1n];for(let i=1;i<=90;i++)factor.push(factor.at(-1)*BigInt(i));
const gcd=(a,b)=>{if(a<0n)a=-a;if(b<0n)b=-b;while(b){let c=a%b;a=b;b=c;}return a;};
const lcm=(a,b)=>a/gcd(a,b)*b;
const dim=72;
const a0=model.A.map(r=>r.map(BigInt));
const b0=model.raw.map(z=>{let n=BigInt(z.num),d=BigInt(z.den);if(n%d)throw Error("bad exact");return n/d;});
const N=model.C.length;
const Ushell=[8,10,12];
const baseTables=[];for(let ss of Ushell){
 let I=new Map();for(let z of model.P.filter(p=>p.s===ss)){let k=z.h.map(Math.abs).join();let t=k+":"+z.i;if(!I.has(t))I.set(t,[z.h.map(Math.abs),z.i,0n]);I.get(t)[2]++;}
 baseTables.push([...I.values()]);
}
const cache=new Map();
function mom(sh,exps){
 let p=[...exps].sort((a,b)=>b-a),key=sh+":"+p.join();if(cache.has(key))return cache.get(key);
 let d=model.U.map(()=>0n),e=0n;
 for(let [h,i,m] of baseTables[sh])d[i]+=m*h.reduce((s,x,j)=>s*BigInt(x)**BigInt(2*p[j]),1n);
 if(sh===0)for(let h of model.roots)e+=h.reduce((s,x,j)=>s*BigInt(x)**BigInt(2*p[j]),1n);
 let r=[d,e];cache.set(key,r);return r;
}
const M3=7;const one=Array.from({length:M3},(_,i)=>i?0n:1n);
const mul=(a,b)=>a.map((_,i)=>a.slice(0,i+1).reduce((z,v,j)=>z+v*b[i-j],0n));
const pow=(a,e)=>{let r=one;for(let i=0;i<e;i++)r=mul(a,r);return r};
const sig=(d,a)=>Array.from({length:d},(_,i)=>i+1).filter(i=>d%i===0).reduce((s,i)=>s+BigInt(i)**BigInt(a),0n);
const ex=a=>one.map((_,i)=>i?(a===3?240n:-504n)*sig(i,a):1n), E4=ex(3),E6=ex(5);
let Del=one.map((_,i)=>i===1?1n:0n);for(let j=1;j<7;j++){let a=one.map((x,i)=>i===j?-1n:x);Del=mul(Del,pow(a,24));}
function nullModular(m){
 const rows=[];let present=[];
 for(let i=0;i<=2;i++){
  let w=2*m-12-12*i,b=(w%4+4)%4===2?1:0,a=(w-6*b)/4;
  if(a>=0){let y=mul(pow(Del,4+i),mul(pow(E4,a),pow(E6,b))).slice(4,7);present.push([i,y]);}
 }
 for(let i=0;i<=2;i++)if(!present.some(x=>x[0]===i)){
  let v=[0n,0n,0n];v[i]=1n;
  for(let [j,y] of present.slice().reverse())if(j<i)v[j]=-(v[0]*y[0]+v[1]*y[1]+v[2]*y[2]);
  rows.push(v);
 }
 return rows;
}
const modularRows=[];
let meta=[];
for(let m=6;m<=24;m++){
 let lambdas=nullModular(m);if(lambdas.length===0)continue;
 let numer=[1n],denom=[1n];
 for(let i=1;i<=m;i++){
  let nu=-numer[i-1]*BigInt((2*m-2*i+2)*(2*m-2*i+1));
  let de=denom[i-1]*BigInt(2*i*(4*m-2*i+dim-2));
  let d=gcd(nu,de);numer.push(nu/d);denom.push(de/d);
 }
 let cden=denom.reduce(lcm,1n),cc=numer.map((z,i)=>z*cden/denom[i]);
 for(let a=m;a>=0;a--)for(let b=Math.min(a,m-a);b>=0;b--){let c=m-a-b;if(c<0||c>b)continue;
  let rr=[],qq=[];
  for(let sh=0;sh<3;sh++){
   let arr=model.U.map(()=>0n), rhs=0n,s=Ushell[sh];
   let pre=cc.map((co,i)=>co*BigInt(16*s)**BigInt(i)*factor[i]*factor[2*m-2*i]);
   for(let ia=0;ia<=a;ia++)for(let ib=0;ib<=b;ib++)for(let ic=0;ic<=c;ic++){
    let i=ia+ib+ic,pp=[a-ia,b-ib,c-ic];
    let w=pre[i]/(factor[ia]*factor[ib]*factor[ic]*pp.reduce((t,x)=>t*factor[2*x],1n));
    let [ro,re]=mom(sh,pp);
    for(let k=0;k<arr.length;k++)arr[k]+=w*ro[k];
    rhs+=w*re;
   }
   rr.push(arr);qq.push(rhs);
  }
  for(let z of lambdas){
   let row=rr[0].map((x,k)=>x*z[0]+rr[1][k]*z[1]+rr[2][k]*z[2]);let rhs=-(z[0]*qq[0]+z[1]*qq[1]+z[2]*qq[2]);let d=row.reduce(gcd,gcd(rhs,0n));if(d===0n)continue;
   modularRows.push([row.map(x=>x/d),rhs/d]);meta.push([m,a,b,c,z.map(String)]);
  }
 }
}
const A=[...a0,...modularRows.map(x=>x[0])],B=[...b0,...modularRows.map(x=>x[1])];
return {A:A.map(r=>r.map(String)),B:B.map(String),meta,modularCount:modularRows.length};
};
const harmD4=function(model){
const factor=[1n];for(let i=1;i<=90;i++)factor.push(factor.at(-1)*BigInt(i));
const gcd=(a,b)=>{if(a<0n)a=-a;if(b<0n)b=-b;while(b){let c=a%b;a=b;b=c;}return a;};
const lcm=(a,b)=>a/gcd(a,b)*b;
const dim=72;
const a0=model.A.map(r=>r.map(BigInt));
const b0=model.raw.map(z=>{let n=BigInt(z.num),d=BigInt(z.den);if(n%d)throw Error("bad exact");return n/d;});
const N=model.C.length;
const Ushell=[8,10,12];
const baseTables=[];for(let ss of Ushell){
 let I=new Map();for(let z of model.P.filter(p=>p.s===ss)){let k=z.h.map(Math.abs).join();let t=k+":"+z.i;if(!I.has(t))I.set(t,[z.h.map(Math.abs),z.i,0n]);I.get(t)[2]++;}
 baseTables.push([...I.values()]);
}
const cache=new Map();
function mom(sh,exps){
 let p=[...exps].sort((a,b)=>b-a),key=sh+":"+p.join();if(cache.has(key))return cache.get(key);
 let d=model.U.map(()=>0n),e=0n;
 for(let [h,i,m] of baseTables[sh])d[i]+=m*h.reduce((s,x,j)=>s*BigInt(x)**BigInt(2*p[j]),1n);
 if(sh===0)for(let h of model.roots)e+=h.reduce((s,x,j)=>s*BigInt(x)**BigInt(2*p[j]),1n);
 let r=[d,e];cache.set(key,r);return r;
}
const M3=7;const one=Array.from({length:M3},(_,i)=>i?0n:1n);
const mul=(a,b)=>a.map((_,i)=>a.slice(0,i+1).reduce((z,v,j)=>z+v*b[i-j],0n));
const pow=(a,e)=>{let r=one;for(let i=0;i<e;i++)r=mul(a,r);return r};
const sig=(d,a)=>Array.from({length:d},(_,i)=>i+1).filter(i=>d%i===0).reduce((s,i)=>s+BigInt(i)**BigInt(a),0n);
const ex=a=>one.map((_,i)=>i?(a===3?240n:-504n)*sig(i,a):1n), E4=ex(3),E6=ex(5);
let Del=one.map((_,i)=>i===1?1n:0n);for(let j=1;j<7;j++){let a=one.map((x,i)=>i===j?-1n:x);Del=mul(Del,pow(a,24));}
function nullModular(m){
 const rows=[];let present=[];
 for(let i=0;i<=2;i++){
  let w=2*m-12-12*i,b=(w%4+4)%4===2?1:0,a=(w-6*b)/4;
  if(a>=0){let y=mul(pow(Del,4+i),mul(pow(E4,a),pow(E6,b))).slice(4,7);present.push([i,y]);}
 }
 for(let i=0;i<=2;i++)if(!present.some(x=>x[0]===i)){
  let v=[0n,0n,0n];v[i]=1n;
  for(let [j,y] of present.slice().reverse())if(j<i)v[j]=-(v[0]*y[0]+v[1]*y[1]+v[2]*y[2]);
  rows.push(v);
 }
 return rows;
}
const modularRows=[];
let meta=[];
for(let m=6;m<=24;m++){
 let lambdas=nullModular(m);if(lambdas.length===0)continue;
 let numer=[1n],denom=[1n];
 for(let i=1;i<=m;i++){
  let nu=-numer[i-1]*BigInt((2*m-2*i+2)*(2*m-2*i+1));
  let de=denom[i-1]*BigInt(2*i*(4*m-2*i+dim-2));
  let d=gcd(nu,de);numer.push(nu/d);denom.push(de/d);
 }
 let cden=denom.reduce(lcm,1n),cc=numer.map((z,i)=>z*cden/denom[i]);
 for(let a=m;a>=0;a--)for(let b=Math.min(a,m-a);b>=0;b--)for(let c=Math.min(b,m-a-b);c>=0;c--){let d=m-a-b-c;if(d<0||d>c)continue;
  let rr=[],qq=[];
  for(let sh=0;sh<3;sh++){
   let arr=model.U.map(()=>0n), rhs=0n,s=Ushell[sh];
   let pre=cc.map((co,i)=>co*BigInt(16*s)**BigInt(i)*factor[i]*factor[2*m-2*i]);
   for(let ia=0;ia<=a;ia++)for(let ib=0;ib<=b;ib++)for(let ic=0;ic<=c;ic++)for(let id=0;id<=d;id++){
    let i=ia+ib+ic+id,pp=[a-ia,b-ib,c-ic,d-id];
    let w=pre[i]/(factor[ia]*factor[ib]*factor[ic]*factor[id]*pp.reduce((t,x)=>t*factor[2*x],1n));
    let [ro,re]=mom(sh,pp);
    for(let k=0;k<arr.length;k++)arr[k]+=w*ro[k];
    rhs+=w*re;
   }
   rr.push(arr);qq.push(rhs);
  }
  for(let z of lambdas){
   let row=rr[0].map((x,k)=>x*z[0]+rr[1][k]*z[1]+rr[2][k]*z[2]);let rhs=-(z[0]*qq[0]+z[1]*qq[1]+z[2]*qq[2]);let d=row.reduce(gcd,gcd(rhs,0n));if(d===0n)continue;
   modularRows.push([row.map(x=>x/d),rhs/d]);meta.push([m,a,b,c,m-a-b-c,z.map(String)]);
  }
 }
}
const A=[...a0,...modularRows.map(x=>x[0])],B=[...b0,...modularRows.map(x=>x[1])];
return {A:A.map(r=>r.map(String)),B:B.map(String),meta,modularCount:modularRows.length};
};
const mod=function(harm){
let A=harm.A.map(r=>r.map(BigInt)),B=harm.B.map(BigInt);const m=A.length,n=A[0].length;
function lap(p){p=BigInt(p);const res=x=>{const y=x%p;return y<0n?y+p:y}; const inv=x=>{let y=p-2n,s=1n;for(;y;y>>=1n,x=x*x%p)if(y&1n)s=s*x%p;return s;};
let W=A.map((r,i)=>[...r.map(res),res(48n*B[i])]);let labels=W.map((_,i)=>i);let ids=[],J=[];
for(let j=0;j<n;j++){let i=W.findIndex((row,i)=>i>=ids.length&&row[j]);if(i<0)continue;let rk=ids.length;[W[i],W[rk]]=[W[rk],W[i]];[labels[i],labels[rk]]=[labels[rk],labels[i]];J.push(labels[rk]);let z=inv(W[rk][j]);for(let k=j;k<=n;k++)W[rk][k]=W[rk][k]*z%p;
for(let v=0;v<W.length;v++){if(v===rk||!W[v][j])continue;let x=W[v][j];for(let k=j;k<=n;k++)W[v][k]=res(W[v][k]-x*W[rk][k]);}ids.push(j);}
let consistent=W.slice(ids.length).every(x=>!x[n]);
const fixed=ids.filter((j,i)=>W[i].slice(0,n).every((x,k)=>!x||ids.includes(k)));
return {p,ids,J,W,consistent,fixed,fv:fixed.map(j=>W[ids.indexOf(j)][n])};}
const out=[1000003].map(lap);
return {out:out.map(x=>({...x,p:String(x.p),W:x.W.map(row=>row.map(String)),fv:x.fv.map(String)}))};
};
const affine=function(harm,rr){
const F=(n,d=1n)=>{if(n===0n)return [0n,1n];if(d<0n){d=-d;n=-n;}let z=G(n,d);return [n/z,d/z]};
const G=(n,d)=>{if(n<0n)n=-n;if(d<0n)d=-d;while(d){let a=n%d;n=d;d=a;}return n};
const plus=([a,b],[c,d])=>F(a*d+b*c,b*d);
const sub=([a,b],[c,d])=>F(a*d-b*c,b*d);
const mul=([a,b],[c,d])=>F(a*c,b*d);
const div=([a,b],[c,d])=>F(a*d,b*c);
let A=harm.A.map(r=>r.map(BigInt)),B=harm.B.map(BigInt);const n=A[0].length;
const rank=rr.ids.length, ix=rr.ids.concat(Array.from({length:n},(_,i)=>i).filter(i=>!rr.ids.includes(i)));
const G0=rr.J.map(j=>ix.map(i=>A[j][i]).concat(B[j]));
let prev=1n,maxl=0;
for(let i=0;i<rank;i++){
 let k=G0.findIndex((z,j)=>j>=i&&z[i]);if(k<0)throw Error("sing "+i);
 [G0[k],G0[i]]=[G0[i],G0[k]];let pivot=G0[i][i];
 for(let j=i+1;j<rank;j++){
  let x=G0[j][i];
  for(let k=i+1;k<=n;k++){
   let t=G0[j][k]*pivot-x*G0[i][k];
   if(t%prev!==0n)throw Error("bareiss");
   G0[j][k]=t/prev;
  }
  G0[j][i]=0n;
 }
 prev=pivot;
}
let gsol=Array.from({length:n},(_,i)=>Array.from({length:n-rank+1},(_,j)=>F(i>=rank&&j===i-rank?1n:0n)));
for(let j=rank-1;j>=0;j--){
 for(let k=0;k<=n-rank;k++){
  let x=k===n-rank?[G0[j][n],1n]:[0n,1n];
  for(let i=j+1;i<n;i++)if(G0[j][i])x=sub(x,mul([G0[j][i],1n],gsol[i][k]));
  gsol[j][k]=div(x,[G0[j][j],1n]);
 }
}
const s=Array.from({length:n},(_,i)=>gsol[ix.indexOf(i)]), zero=[0n,1n];
for(let i=0;i<A.length;i++){
  for(let k=0;k<=n-rank;k++){
   let c=s.reduce((z,x,j)=>A[i][j]?plus(z,mul(x[k],[A[i][j],1n])):z,zero);
   if(c[1]!==1n||c[0]!== (k===n-rank?B[i]:0n))throw Error("bad "+JSON.stringify([i,k,c.map(String)]));
 }}
return {s:s.map(r=>r.map(x=>x.map(String))),free:ix.slice(rank)};
};
const exactB=function(e){
const p=(a,b)=>a+b;
const G=(n,d)=>{if(n<0n)n=-n;if(d<0n)d=-d;while(d){let a=n%d;n=d;d=a;}return n};
const F=(n,d=1n)=>{if(n===0n)return [0n,1n];if(d===0n)throw Error("0-den");if(d<0n){d=-d;n=-n;}let z=G(n,d);return [n/z,d/z]};
const plus=([a,b],[c,d])=>F(a*d+b*c,b*d);
const sub=([a,b],[c,d])=>F(a*d-b*c,b*d);
const mul=([a,b],[c,d])=>F(a*c,b*d);
const div=([a,b],[c,d])=>F(a*d,b*c);
const S=e.s.map(r=>r.map(x=>x.map(BigInt))),n=e.free.length;
function exactValue(j,sgn,act){
 const a=S[j].map(([x,d])=>[BigInt(sgn)*x,d]); // I want proof n_j for lower, -n_j for upper.
 let M=Array.from({length:n},(_,i)=>act.map(r=>S[r][i]).concat([a[i]]));
 let k=0,J=[];for(let i=0;i<act.length;i++){
  let j=M.findIndex((r,p)=>p>=k&&r[i][0]);if(j<0)continue;
  [M[k],M[j]]=[M[j],M[k]];const x=M[k][i];M[k]=M[k].map(z=>div(z,x));
  for(let j=0;j<M.length;j++)if(j!==k&&M[j][i][0]){
   let z=M[j][i];M[j]=M[j].map((x,l)=>sub(x,mul(z,M[k][l])));
  }
  J.push(i);k++;
 }
 if(!M.slice(k).every(row=>!row.at(-1)[0]))return null;
 let q=act.map((_,i)=>J.includes(i)?M[J.indexOf(i)].at(-1):[0n,1n]);
 if(q.some(x=>x[0]<0n))return null;
 let dif=a.map((v,i)=>sub(v,act.reduce((s,j,k)=>plus(s,mul(q[k],S[j][i])),[0n,1n])));
 if(dif.slice(0,n).some(v=>v[0]))return null;
 let bound=dif[n];if(sgn===-1)bound=[-bound[0],bound[1]];
 return {i:j,sgn,bound:bound.map(String),expr:act.map((k,i)=>[k,q[i].map(String)]).filter(x=>x[1][0]!=="0")};
}
return exactValue;
};
const heavy=function(m){
const C=m.C, ids=C.map((v,i)=>i).filter(i=>C[i].r===124); 
function dist(a,b){let norm=0,k=0,qmax=0;for(let i=0;i<4;i++){let d=a.h[i]-b.h[i],r=Math.floor((d+4)/8),z=d-r*8,az=Math.abs(z);k+=r;norm+=z*z;qmax=Math.max(qmax,az);}return norm+((k&1)?64-16*qmax:0);}
const E=ids.map(i=>ids.map(j=>dist(C[i],C[j])));
const deg=E.map(r=>r.filter(x=>x>0&&x<=16).length);
const choices=ids.map((i,p)=>({i,p})).filter(({i,p})=>!ids.slice(0,p).some(j=>C[j].i===C[i].i));

function perm(a){if(!a.length)return [[]];let z=[];for(let i=0;i<a.length;i++)for(let b of perm(a.filter((_,j)=>i!==j)))z.push([a[i],...b]);return z;}
function feasible(S){
 const k=S.length,lookup=S.map((a,i)=>S.map(b=>E[a][b]));
 if(k<3)return true;
 for(let mask=0;mask<(1<<(k-1));mask++){
  let a=[],b=[]; for(let i=0;i<k;i++)((mask>>i)&1?b:a).push(i);
  if(a.some(i=>a.some(j=>lookup[i][j]>4))||b.some(i=>b.some(j=>lookup[i][j]>4)))continue;
  for(let A of perm(a))for(let B of perm(b))if(A.every((i,r)=>B.every((j,s)=>lookup[i][j] <= (r+s===0?16:r+s===1?8:4))))return true;
 }
 return false;
}

let maxima=[],Bn=0;
for(let {i,p} of choices){
 let Q=ids.map((_,j)=>j).filter(j=>j!==p&&E[p][j]>0&&E[p][j]<=16);
 let best=0,wits=[],stats=Array(10).fill(0),xbest=0,Xw=[];
 function ac(P,S){
  Bn++;
  const k=S.length;
  if(k>best){best=k;wits=S;}
  stats[k]++;
  if(k>=8)return;
  // Necessary exact sign/endpoint constraints; delta 16 means both extreme at opposite heights.
  let adj=S.map((i,a)=>S.map((j,b)=>E[i][j]>=8?b:-1).filter(j=>j>=0));
  let color=Array(k).fill(null),comp=Array(k).fill(null),t=0;
  for(let i=0;i<k;i++)if(color[i]===null){let todo=[i];color[i]=0;comp[i]=t++;while(todo.length){let j=todo.pop();for(let b of adj[j])if(color[b]===null){color[b]=1-color[j];comp[b]=comp[j];todo.push(b);}else if(color[b]===color[j])return;}}
  let ext=S.map(i=>S.some(j=>E[i][j]===16));
  for(let a=0;a<k;a++)for(let b=0;b<a;b++)if(ext[a]&&ext[b]&&comp[a]===comp[b]&&color[a]===color[b])return;
  // Across components, two endpoint-mandatory vertices must be opposite regardless of closeness.
  let exi=ext.filter(Boolean).length;if(exi>2)return;
  if(!feasible(S))return;
  if(k>xbest){xbest=k;Xw=S;}
  for(let j=0;j<P.length;j++){
   let u=P[j],can=P.slice(j+1).filter(v=>E[u][v]<=16);
   ac(can,[...S,u]);
  }
 }
 ac(Q,[p]);maxima.push({rep:C[i],deg:Q.length,stats,geomOnly:best,signBound:xbest,X:Xw.map(j=>C[ids[j]])});
}
return {n:ids.length,maxDeg:Math.max(...deg),Bn,maxima};
};
const relaxed=function(m,E0){
const C=m.C, ids=C.map((v,i)=>i).filter(i=>C[i].r===124); 
function dist(a,b){let norm=0,k=0,qmax=0;for(let i=0;i<4;i++){let d=a.h[i]-b.h[i],r=Math.floor((d+4)/8),z=d-r*8,az=Math.abs(z);k+=r;norm+=z*z;qmax=Math.max(qmax,az);}return norm+((k&1)?64-16*qmax:0);}
const E=E0;
const deg=E.map(r=>r.filter(x=>x>0&&x<=16).length);
const choices=ids.map((i,p)=>({i,p})).filter(({i,p})=>!ids.slice(0,p).some(j=>C[j].i===C[i].i));

function perm(a){if(!a.length)return [[]];let z=[];for(let i=0;i<a.length;i++)for(let b of perm(a.filter((_,j)=>i!==j)))z.push([a[i],...b]);return z;}
function feasible(S){
 const k=S.length,lookup=S.map((a,i)=>S.map(b=>E[a][b]));
 if(k<3)return true;
 for(let mask=0;mask<(1<<(k-1));mask++){
  let a=[],b=[]; for(let i=0;i<k;i++)((mask>>i)&1?b:a).push(i);
  if(a.some(i=>a.some(j=>lookup[i][j]>4))||b.some(i=>b.some(j=>lookup[i][j]>4)))continue;
  const depth=L=>{let d=[];for(let i=0;i<L.length;i++){let r=d.at(-1)??0;for(let j=0;j<i;j++)if(lookup[L[i]][L[j]]===4)r=Math.max(r,1+d[j]);d.push(r);}return d;};
 for(let A of perm(a)){let da=depth(A);for(let B of perm(b)){let db=depth(B);if(A.every((i,r)=>B.every((j,s)=>lookup[i][j] <= (da[r]+db[s]===0?16:da[r]+db[s]===1?8:4))))return true;}}
 }
 return false;
}

let maxima=[],Bn=0;
for(let {i,p} of choices){
 let Q=ids.map((_,j)=>j).filter(j=>j!==p&&E[p][j]<=16);
 let best=0,wits=[],stats=Array(10).fill(0),xbest=0,Xw=[];
 function ac(P,S){
  Bn++;
  const k=S.length;
  if(k>best){best=k;wits=S;}
  stats[k]++;
  if(k>=8)return;
  // Necessary exact sign/endpoint constraints; delta 16 means both extreme at opposite heights.
  let adj=S.map((i,a)=>S.map((j,b)=>E[i][j]>=8?b:-1).filter(j=>j>=0));
  let color=Array(k).fill(null),comp=Array(k).fill(null),t=0;
  for(let i=0;i<k;i++)if(color[i]===null){let todo=[i];color[i]=0;comp[i]=t++;while(todo.length){let j=todo.pop();for(let b of adj[j])if(color[b]===null){color[b]=1-color[j];comp[b]=comp[j];todo.push(b);}else if(color[b]===color[j])return;}}
  let ext=S.map(i=>S.some(j=>E[i][j]===16));
  for(let a=0;a<k;a++)for(let b=0;b<a;b++)if(ext[a]&&ext[b]&&E[S[a]][S[b]]>0&&comp[a]===comp[b]&&color[a]===color[b])return;
  // Across components, two endpoint-mandatory vertices must be opposite regardless of closeness.
  
  if(!feasible(S))return;
  if(k>xbest){xbest=k;Xw=S;}
  for(let j=0;j<P.length;j++){
   let u=P[j],can=P.slice(j+1).filter(v=>E[u][v]<=16);
   ac(can,[...S,u]);
  }
 }
 ac(Q,[p]);maxima.push({rep:C[i],deg:Q.length,stats,geomOnly:best,signBound:xbest,X:Xw.map(j=>C[ids[j]])});
}
return {n:ids.length,maxDeg:Math.max(...deg),Bn,maxima};
};

function need(x,w){if(!x)throw Error(w);}
{
 const pl=(a,b)=>a**4n-126n*a*a*b*b+b**4n;
 need(pl(89n,1000n)>0n&&pl(90n,1000n)<0n,"U enclosure");
 const pr=(a,b)=>a*a-4n*a*b+b*b;
 need(pr(267n,1000n)>0n&&pr(268n,1000n)<0n,"R enclosure");
 const A=[[89n,1000n],[90n,1000n]],B=[[267n,1000n],[268n,1000n]];
 const mul=(x,y)=>[x[0]*y[0],x[1]*y[1]];
 function cmp([a,b],[c,d],threshold){
  const nn=a*a*d*d+c*c*b*b+4n*a*b*c*d, dd=(b*b-a*a)*(d*d-c*c);
  return 248n*nn-BigInt(threshold-4)*dd;
 }
 need(cmp(A[0],mul(A[0],B[0]),8)>0n,"mixed rank lower");
 need(cmp(A[1],mul(A[1],B[1]),12)<0n,"mixed rank upper");
 need(cmp(A[1],mul(A[1],mul(B[1],B[1])),8)<0n,"rank 0/2 upper");
 need(cmp(mul(A[1],B[1]),mul(A[1],B[1]),8)<0n,"rank 1/1 upper");
 // Same sign with raw squared distance ≥ 8: c≥16/31. The four
 // rectangle corners of the separately concave dot all exceed one half.
 need(2*16>31&&4*8>31&&2*(16+1)>32,"same sign corners");
}

const rid=[0,3,4,16],G=SOURCE.catalogue_gram;
need(G.length===72&&G.every(g=>g.length===72),"catalogue shape");
need(JSON.stringify(rid.map(i=>rid.map(j=>G[i][j])))==="[[8,-4,-4,-4],[-4,8,0,0],[-4,0,8,0],[-4,0,0,8]]","explicit D4 Gram");
const a=makeA3(), b=makeD4(a.theta);
need(a.theta.slice(0,7).join()==="1,0,0,0,6218175600,15281788354560,9026867482214400","theta");
function solve(m,hs,norm,rank){
 need(m.U.length===norm,"class catalogue");
 const h=hs(m), col=m.U.map((_,i)=>i).sort((i,j)=>m.U[j].s0-m.U[i].s0||m.U[j].r-m.U[i].r);
 const dat={A:h.A.map(r=>col.map(i=>r[i])),B:h.B};
 const rr=mod(dat).out[0];
 need(rr.consistent&&rr.ids.length===rank,"rank/consistency");
 return affine(h,{...rr,ids:rr.ids.map(i=>col[i])});
}
const ea=solve(a,harmA3,66,61),eb=solve(b,harmD4,48,44);
need(ea.free.join()==="0,1,2,5,8"&&eb.free.join()==="0,1,2,3","free population axes");
for(const [i,sgn,act,want] of [[44,1,[3,0,1,5,8],"627133872"],[44,-1,[3,9,22,1,5],"642626560"]]){
 const z=exactB(ea)(i,sgn,act);
 need(z&&z.bound.join()===want+",1","A3 deep bound");
}
const cuts=[{"i":0,"sgn":1,"bound":["0","1"],"expr":[[0,["1","1"]]]},{"i":0,"sgn":-1,"bound":["1160","1"],"expr":[[4,["3","20"]],[16,["1","80"]],[1,["16","5"]],[3,["8","5"]]]},{"i":1,"sgn":1,"bound":["0","1"],"expr":[[1,["1","1"]]]},{"i":1,"sgn":-1,"bound":["337","1"],"expr":[[4,["3","32"]],[0,["1","8"]],[2,["1","8"]],[3,["1","2"]]]},{"i":2,"sgn":1,"bound":["0","1"],"expr":[[2,["1","1"]]]},{"i":2,"sgn":-1,"bound":["2696","1"],"expr":[[4,["3","4"]],[0,["1","1"]],[1,["8","1"]],[3,["4","1"]]]},{"i":3,"sgn":1,"bound":["0","1"],"expr":[[3,["1","1"]]]},{"i":3,"sgn":-1,"bound":["674","1"],"expr":[[4,["3","16"]],[0,["1","4"]],[1,["2","1"]],[2,["1","4"]]]},{"i":4,"sgn":1,"bound":["0","1"],"expr":[[4,["1","1"]]]},{"i":4,"sgn":-1,"bound":["10784","3"],"expr":[[0,["4","3"]],[1,["32","3"]],[2,["4","3"]],[3,["16","3"]]]},{"i":5,"sgn":1,"bound":["0","1"],"expr":[[3,["20","1"]]]},{"i":5,"sgn":-1,"bound":["13480","1"],"expr":[[4,["15","4"]],[0,["5","1"]],[1,["40","1"]],[2,["5","1"]]]},{"i":6,"sgn":1,"bound":["0","1"],"expr":[[4,["2","1"]],[1,["22","1"]],[3,["10","1"]]]},{"i":6,"sgn":-1,"bound":["7414","1"],"expr":[[4,["1","16"]],[0,["11","4"]],[2,["11","4"]],[3,["1","1"]]]},{"i":7,"sgn":1,"bound":["87040","1"],"expr":[[16,["3","8"]],[1,["64","1"]],[2,["2","1"]],[3,["23","1"]]]},{"i":7,"sgn":-1,"bound":["127232","1"],"expr":[[4,["6","1"]],[0,["32","1"]],[1,["48","1"]],[3,["33","1"]]]},{"i":8,"sgn":1,"bound":["19952","1"],"expr":[[0,["34","1"]],[1,["32","1"]],[2,["28","1"]],[3,["16","1"]]]},{"i":8,"sgn":-1,"bound":["102400","1"],"expr":[[4,["219","10"]],[16,["3","40"]],[1,["1056","5"]],[3,["528","5"]]]},{"i":9,"sgn":1,"bound":["380896","1"],"expr":[[0,["196","1"]],[1,["32","1"]],[2,["4","1"]],[3,["16","1"]]]},{"i":9,"sgn":-1,"bound":["614400","1"],"expr":[[4,["159","5"]],[16,["12","5"]],[1,["3072","5"]],[3,["1536","5"]]]},{"i":10,"sgn":1,"bound":["1146880","1"],"expr":[[4,["176","5"]],[16,["8","5"]],[1,["2548","5"]],[3,["1364","5"]]]},{"i":10,"sgn":-1,"bound":["4058624","3"],"expr":[[0,["448","3"]],[1,["212","3"]],[2,["64","3"]],[3,["52","3"]]]},{"i":11,"sgn":1,"bound":["4295200","1"],"expr":[[4,["51","1"]],[0,["1020","1"]],[1,["736","1"]],[3,["240","1"]]]},{"i":11,"sgn":-1,"bound":["5478400","1"],"expr":[[4,["102","1"]],[16,["51","4"]],[1,["2528","1"]],[3,["1392","1"]]]},{"i":12,"sgn":1,"bound":["13516800","1"],"expr":[[4,["348","5"]],[16,["44","5"]],[1,["8544","5"]],[3,["4827","5"]]]},{"i":12,"sgn":-1,"bound":["14333440","1"],"expr":[[4,["36","1"]],[0,["704","1"]],[1,["544","1"]],[3,["161","1"]]]},{"i":13,"sgn":1,"bound":["43485184","1"],"expr":[[4,["6","1"]],[0,["3840","1"]],[1,["1536","1"]]]},{"i":13,"sgn":-1,"bound":["47939584","1"],"expr":[[4,["570","1"]],[16,["48","1"]],[1,["10752","1"]],[3,["6144","1"]]]},{"i":14,"sgn":1,"bound":["125501440","1"],"expr":[[4,["2088","5"]],[16,["144","5"]],[1,["35064","5"]],[3,["20172","5"]]]},{"i":14,"sgn":-1,"bound":["128432896","1"],"expr":[[0,["2400","1"]],[1,["1128","1"]],[2,["96","1"]],[3,["36","1"]]]},{"i":15,"sgn":1,"bound":["365817816","1"],"expr":[[0,["13176","1"]],[1,["6912","1"]],[2,["864","1"]],[3,["384","1"]]]},{"i":15,"sgn":-1,"bound":["382429080","1"],"expr":[[4,["12474","5"]],[16,["1539","10"]],[1,["196992","5"]],[3,["113856","5"]]]},{"i":16,"sgn":1,"bound":["0","1"],"expr":[[16,["1","1"]]]},{"i":16,"sgn":-1,"bound":["92800","1"],"expr":[[4,["12","1"]],[0,["80","1"]],[1,["256","1"]],[3,["128","1"]]]},{"i":17,"sgn":1,"bound":["1146880","1"],"expr":[[4,["156","5"]],[16,["8","5"]],[1,["3363","5"]],[3,["1364","5"]]]},{"i":17,"sgn":-1,"bound":["1383991","1"],"expr":[[4,["405","32"]],[0,["1287","8"]],[2,["263","8"]],[3,["127","2"]]]},{"i":18,"sgn":1,"bound":["4085120","1"],"expr":[[4,["102","1"]],[0,["1104","1"]],[1,["896","1"]],[3,["576","1"]]]},{"i":18,"sgn":-1,"bound":["5365760","1"],"expr":[[4,["318","5"]],[16,["69","5"]],[1,["13184","5"]],[3,["5952","5"]]]},{"i":19,"sgn":1,"bound":["13516800","1"],"expr":[[4,["348","5"]],[16,["44","5"]],[1,["8544","5"]],[3,["4162","5"]]]},{"i":19,"sgn":-1,"bound":["14333440","1"],"expr":[[4,["36","1"]],[0,["704","1"]],[1,["544","1"]],[3,["294","1"]]]},{"i":20,"sgn":1,"bound":["129765376","3"],"expr":[[0,["11776","3"]],[1,["512","3"]],[2,["256","3"]],[3,["1024","3"]]]},{"i":20,"sgn":-1,"bound":["47939584","1"],"expr":[[4,["640","1"]],[16,["48","1"]],[1,["12800","1"]],[3,["6144","1"]]]},{"i":21,"sgn":1,"bound":["125501440","1"],"expr":[[4,["2088","5"]],[16,["144","5"]],[1,["35064","5"]],[3,["14352","5"]]]},{"i":21,"sgn":-1,"bound":["128432896","1"],"expr":[[0,["2400","1"]],[1,["1128","1"]],[2,["96","1"]],[3,["1200","1"]]]},{"i":22,"sgn":1,"bound":["125501440","1"],"expr":[[4,["1438","5"]],[16,["144","5"]],[1,["35394","5"]],[3,["17262","5"]]]},{"i":22,"sgn":-1,"bound":["128174080","1"],"expr":[[4,["58","1"]],[0,["2304","1"]],[1,["294","1"]],[3,["234","1"]]]},{"i":24,"sgn":1,"bound":["362622720","1"],"expr":[[4,["768","1"]],[0,["11808","1"]],[1,["6144","1"]],[3,["3072","1"]]]},{"i":24,"sgn":-1,"bound":["376320000","1"],"expr":[[4,["5016","5"]],[16,["738","5"]],[1,["158208","5"]],[3,["79104","5"]]]},{"i":32,"sgn":1,"bound":["362802880","1"],"expr":[[4,["216","1"]],[0,["15400","1"]],[1,["8192","1"]],[3,["5120","1"]]]},{"i":32,"sgn":-1,"bound":["380666880","1"],"expr":[[4,["2094","1"]],[16,["385","2"]],[1,["41088","1"]],[3,["19520","1"]]]}];
for(const cut of cuts){
 const z=exactB(eb)(cut.i,cut.sgn,cut.expr.map(x=>x[0]));
 need(z&&z.bound.join()===cut.bound.join(),"D4 affine cut");
}
const plain=heavy(b), C=b.C.filter(v=>v.r===124);
const E=C.map(u=>C.map(v=>{const q=b.minQ(u.h.map((x,i)=>x-v.h[i]));return q<36?q:q===64?0:q-32;}));
const loose=relaxed(b,E);
need(C.length===192&&plain.maxima.map(x=>x.signBound).join()==="5,2,4","strict heavy cap");
need(loose.maxima.map(x=>x.signBound).join()==="5,5,5","fictional shorter-lattice improvements");
const ret = {
 A3:{classes:66,rank:61,deep_hole_norm12:[627133872,642626560]},
 D4:{classes:48,rank:44,projection_cosets:b.C.filter(x=>x.s===8).length, heavy:C.length,heavy_caps_all:[5,2,4],heavy_caps_hypothetical:[5,5,5],weights:[
   {source_norm:10,radius2:8,bounds:[362622720,376320000]},
   {source_norm:12,radius2:8,bounds:[362802880,380666880]},
   {source_norm:10,radius2:31/4,bounds:[125501440,128432896]}]},
 scope:"A3/D4 populations bound orientation averages, not all actual fibres. Heavy caps refer to the universal wrapped shortest-distance interface and its stated artificial relaxation; not an unrestricted upper bound."
};
if(typeof report==="function")report(ret);else console.log(JSON.stringify(ret,null,2));

