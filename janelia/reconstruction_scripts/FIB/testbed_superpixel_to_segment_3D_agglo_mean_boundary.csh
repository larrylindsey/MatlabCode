# ./superpixel_to_segment_3D_agglo_mean_boundary.csh "rat_cortex.FIB.013009" 46 75 mf5 ".ws.T60" ".ld.T60_L1000" 70 80 90 100 110 120 130 140 150 160

# ./blend_stack_segment.csh "rat_cortex.FIB.013009" 46 75 mf5 ".ws.T60" "ladder" ".ld.T60_L1000" ~/temp/fib.test_ladder.46.75.tif

./blend_stack_segment.csh "rat_cortex.FIB.013009" 46 75 mf5 ".ws.T60" "agglo_mean_boundary" ".ld.T60_L1000.amb.T140" ~/temp/fib.test_amb.T140.46.75.tif 
