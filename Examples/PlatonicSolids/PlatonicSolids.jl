push!(LOAD_PATH, "Modules")
using Colors
using JuliRay

φ=MathConstants.φ

## Tetrahedron
vertices=[[1.0,1.0,1.0],[1.0,-1.0,-1.0],[-1.0,1.0,-1.0],[-1.0,-1.0,1.0]]
edges=[[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
faces=[[2,3,4],[1,3,4],[1,2,4],[1,2,3]]

V=rgbColor(csgUnion((p->Sphere(p,0.05)).(vertices)),RGB(0.1,0.1,0.1))
E=rgbColor(csgUnion([Cylinder(vertices[I...],vertices[J...],0.025) for (I,J) in edges]),RGB(0.2,0.2,0.2))
F=rgbftColor(csgUnion([Polygon((i->vertices[i]).(I)) for I in faces]),RGB(1,0,0),FT(0.2,0.2))
render(csgUnion(V,E,F),camera=LngLatCamera(lng=π/5,lat=π/5,pers=0.2,zoom=0.25),name="PlatonicSolid_4")

## Cube
vertices=[[[i,j,k] for i ∈ [1,-1], j ∈ [1,-1], k ∈ [1,-1]]...]
edges=[[[1+i+j,2+i+j] for i ∈ [0,2], j ∈ [0,4]]...]
edges=vcat(edges,[[[1+i+j,3+i+j] for i ∈ [0,1], j ∈ [0,4]]...])
edges=vcat(edges,[[[1+i+j,5+i+j] for i ∈ [0,1], j ∈ [0,2]]...])
faces=[[1,2,4,3],[1,2,4,3].+4,[1,3,7,5],[1,3,7,5].+1,[1,5,6,2],[1,5,6,2].+2]

V=rgbColor(csgUnion((p->Sphere(p,0.05)).(vertices)),RGB(0.1,0.1,0.1))
E=rgbColor(csgUnion([Cylinder(vertices[I...],vertices[J...],0.025) for (I,J) in edges]),RGB(0.2,0.2,0.2))
F=rgbftColor(csgUnion([Polygon((i->vertices[i]).(I)) for I in faces]),RGB(1,0,0),FT(0.2,0.2))
render(csgUnion(V,E,F),camera=LngLatCamera(lng=π/5,lat=π/5,pers=0.2,zoom=0.25),name="PlatonicSolid_6")

## Octahedron
vertices=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]*sqrt(3)
edges=[[1,3],[3,2],[2,4],[4,1],[3,5],[5,4],[4,6],[6,3],[5,1],[1,6],[6,2],[2,5]]
faces=[[1,3,5],[2,3,5],[2,4,5],[1,4,5],[1,3,6],[2,3,6],[2,4,6],[1,4,6]]

V=rgbColor(csgUnion((p->Sphere(p,0.05)).(vertices)),RGB(0.1,0.1,0.1))
E=rgbColor(csgUnion([Cylinder(vertices[I...],vertices[J...],0.025) for (I,J) in edges]),RGB(0.2,0.2,0.2))
F=rgbftColor(csgUnion([Polygon((i->vertices[i]).(I)) for I in faces]),RGB(1,0,0),FT(0.2,0.2))
render(csgUnion(V,E,F),camera=LngLatCamera(lng=π/5,lat=π/5,pers=0.2,zoom=0.25),name="PlatonicSolid_8")

## Dodecahedron
vertices0=[[[i,j,k] for i ∈ [1,-1], j ∈ [1,-1], k ∈ [1,-1]]...]
vertices1=[[[0,i/φ,j*φ] for i ∈ [1,-1], j ∈ [1,-1]]...]
vertices2=[[[j*φ,0,i/φ] for i ∈ [1,-1], j ∈ [1,-1]]...]
vertices3=[[[i/φ,j*φ,0] for i ∈ [1,-1], j ∈ [1,-1]]...]
vertices=vcat(vertices0,vertices1,vertices2,vertices3)
edges1=[[1,9],[2,9],[3,10],[4,10],[9,10]]
edges2=[[5,11],[6,11],[7,12],[8,12],[11,12]]
edges3=[[1,13],[3,13],[5,14],[7,14],[13,14]]
edges4=[[2,15],[4,15],[6,16],[8,16],[15,16]]
edges5=[[1,17],[5,17],[2,18],[6,18],[17,18]]
edges6=[[3,19],[7,19],[4,20],[8,20],[19,20]]
edges=vcat(edges1,edges2,edges3,edges4,edges5,edges6)
faces=[[1,9,10,3,13],[2,9,10,4,15],[5,11,12,7,14],[6,11,12,8,16],[1,13,14,5,17],[3,13,14,7,19],[2,15,16,6,18],[4,15,16,8,20],[1,17,18,2,9],[5,17,18,6,11],[3,19,20,4,10],[7,19,20,8,12]]

V=rgbColor(csgUnion((p->Sphere(p,0.05)).(vertices)),RGB(0.1,0.1,0.1))
E=rgbColor(csgUnion([Cylinder(vertices[I...],vertices[J...],0.025) for (I,J) in edges]),RGB(0.2,0.2,0.2))
F=rgbftColor(csgUnion([Polygon((i->vertices[i]).(I)) for I in faces]),RGB(1,0,0),FT(0.2,0.2))
render(csgUnion(V,E,F),camera=LngLatCamera(lng=π/5,lat=π/5,pers=0.2,zoom=0.25),name="PlatonicSolid_12")

## Icosahedron
vertices1=[[[0,i,j*φ] for i ∈ [1,-1], j ∈ [1,-1]]...]
vertices2=[[[j*φ,0,i] for i ∈ [1,-1], j ∈ [1,-1]]...]
vertices3=[[[i,j*φ,0] for i ∈ [1,-1], j ∈ [1,-1]]...]
vertices=vcat(vertices1,vertices2,vertices3)*sqrt(3)/√(φ+2)
edges=[[2i-1,2i] for i ∈ 1:6]
edges=[edges...,[1,5],[5,9],[9,1],[2,5],[5,11],[11,2],[3,6],[6,9],[9,3],[4,6],[6,11],[11,4],[1,7],[7,10],[10,1],[2,7],[7,12],[12,2],[3,8],[8,10],[10,3],[4,8],[8,12],[12,4]]
faces=[[1,5,9],[2,5,11],[3,6,9],[4,6,11],[1,7,10],[2,7,12],[3,8,10],[4,8,12],[1,2,5],[1,2,7],[3,4,6],[3,4,8],[5,6,9],[5,6,11],[7,8,10],[7,8,12],[9,10,1],[9,10,3],[11,12,2],[11,12,4]]

V=rgbColor(csgUnion((p->Sphere(p,0.05)).(vertices)),RGB(0.1,0.1,0.1))
E=rgbColor(csgUnion([Cylinder(vertices[I...],vertices[J...],0.025) for (I,J) in edges]),RGB(0.2,0.2,0.2))
F=rgbftColor(csgUnion([Polygon((i->vertices[i]).(I)) for I in faces]),RGB(1,0,0),FT(0.2,0.2))
render(csgUnion(V,E,F),camera=LngLatCamera(lng=π/5,lat=π/5,pers=0.2,zoom=0.25),name="PlatonicSolid_20")
