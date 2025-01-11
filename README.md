## 代码运行方法
 - 直接运行driver_quradrangle.m和driver_triangle.m两个文件，可分别进行四边形网格和三角形网格的有限元分析  
 - 修改quarter-plate-with-hole-quad.geo和quarter-plate-with-hole-tria.geo中Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 后的数字，可改变网格密度，在Gmesh中导出msh文件，文件版本选择version 2 ASCII，只勾选save all element，覆盖同名文件后即可再次运行。  
 - driver_quradrangle.m的第198、199和204、205行可通过添加或取消注释边界条件修改右边和上边施加的应力边界条件  
 - Verification_quad.m为组装解的测试，但无法给出正确的结果
