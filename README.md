Features of the code:-
- Allows the user to interactively select an image (Of type - jpg, png, bmp,tif) from a folder.
- The selected image is loaded and processed for shape matching using morphological operations and filtering.
- The code extracts red pixels from the original image based on specified HSV criteria.
- It applies edge detection to obtain a binary edge map.
- The binary edge map is further processed to create a filled polygon shape.
- The code then identifies the three largest objects in the binary mask.
- A region of interest is extracted from the original image based on the identified largest objects to improve efficiency.
- Searches for the best-fit polygon shape by varying position and scale.
- The results are displayed using a figure with subplots showing various image processing stages and the best-fit shape overlaid on the processed image.
