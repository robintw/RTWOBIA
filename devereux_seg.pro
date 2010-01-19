PRO DEVEREUX_SEG, t
  ENVI_SELECT, fid=fid, dims=dims, pos=pos
  
  image = ENVI_GET_DATA(fid=fid, dims=dims, pos=pos)
  
  result = CREATE_EDGE_MAP(fid, dims, pos, 15.0)
  
  ENVI_ENTER_DATA, result
  
  ;edge_map = ROBERTS(image)
  
  ; Change this bit later to do proper thresholding
  ;edge_map[WHERE(edge_map LE 15)] = 0
  ;edge_map[WHERE(edge_map GT 15)] = 1
  
  ;local_max = intarr(dims[2], dims[4])
  ;indices_of_max = PICKMIN(edge_map, window=81, /max)
  
  
  
  ;local_max[indices_of_max] = 1
  
  ;ENVI_ENTER_DATA, local_max
  
END

FUNCTION SELECT_LOCAL_MAXIMA, array
  s = SIZE(array, /DIMENSIONS)
  
  output = intarr(s[0], s[1])
  
  FOR i = 0, s[0] DO BEGIN
    FOR j = 0, s[1] DO BEGIN
      
    ENDFOR
  ENDFOR
END

FUNCTION DO_CONVOL, x, y, nb, array
  kernel = intarr(3, 3)
  kernel[x, y] = 1
  
  ;result = array
  dims = SIZE(array, /DIMENSIONS)
  result = fltarr(dims[0], dims[1], dims[2], /NOZERO)
  
  
  
  
  FOR i = 0, nb-1 DO BEGIN
    result[*, *, i] = CONVOL(float(array[*, *, i]), kernel, /CENTER, /EDGE_TRUNCATE)
  ENDFOR
  
  return, result
END

FUNCTION CREATE_EDGE_MAP, fid, dims, pos, percentage
  COMPILE_OPT STRICTARR
  
  proj = ENVI_GET_PROJECTION(fid=fid, pixel_size=pixel_size)
  pixel_size = pixel_size[0]
  
  ENVI_FILE_QUERY, fid, nb=nb
 
  ; Get the whole image cube
  
  tile_id = ENVI_INIT_TILE(fid, pos, interleave=0, num_tiles=num_tiles, xs=dims[1], xe=dims[2], ys=dims[3], ye=dims[4])
  
  WholeImage = intarr(dims[2] + 1, dims[4] + 1, nb)
  
  FOR i = 0, num_tiles-1 DO BEGIN
    data = ENVI_GET_TILE(tile_id, i, band_index=band_index)
    WholeImage[*, *, i] = ENVI_GET_TILE(tile_id, i, band_index=band_index)
  ENDFOR
  
  ENVI_TILE_DONE, tile_id

  Z1 = DO_CONVOL(0, 0, nb, WholeImage)
  Z2 = DO_CONVOL(1, 0, nb, WholeImage)
  Z3 = DO_CONVOL(2, 0, nb, WholeImage)
  Z4 = DO_CONVOL(0, 1, nb, WholeImage)
  Z5 = DO_CONVOL(1, 1, nb, WholeImage)
  Z6 = DO_CONVOL(2, 1, nb, WholeImage)
  Z7 = DO_CONVOL(0, 2, nb, WholeImage)
  Z8 = DO_CONVOL(1, 2, nb, WholeImage)
  Z9 = DO_CONVOL(2, 2, nb, WholeImage)
  
  r = fltarr(dims[2] + 1, dims[4] + 1)
  t = fltarr(dims[2] + 1, dims[4] + 1)
  
  
  ; Calculate the partial derivatives
  r = ((Z1 + Z3 + Z4 + Z6 + Z7 + Z9) - (2 * (Z2 + Z5 + Z8))) / (3 * pixel_size^2)
  
  t = ((Z1 + Z2 + Z3 + Z7 + Z8 + Z9) - (2 * (Z4 + Z5 + Z6))) / (3 * pixel_size^2)
  
  t_plus_r = t + r
  
  sums_of_bands = TOTAL(t_plus_r, nb)
  
  ; Calculate the multispectral slope, and replace all NaNs with zero
  ms_slope = sqrt(sums_of_bands)
  
  nan_indices = WHERE(FINITE(ms_slope) EQ 0)
  ms_slope[nan_indices] = 0
  
  ; Calculate the threshold for the slope image
  sorted_image_indices = REVERSE(SORT(ms_slope))
  threshold =  ms_slope[sorted_image_indices[float(percentage)/100 * N_ELEMENTS(ms_slope)]]
  
  print, threshold
  
  ; Create the edge image by thresholding
  edge_image = intarr(dims[2] + 1, dims[4] + 1)
  
  edge_image[WHERE(ms_slope GT threshold)] = 1
  
  return, edge_image  
END

FUNCTION INEQUALITY_FUNCTION, array, t
  result = abs(array - mean(array))^2 / t
  
  return, total(result)
END