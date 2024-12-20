## Huffman Coding Compression

### Problem Description:
Design and implement a file compression program using the Huffman Coding algorithm for a given PGM image file. 

#### Requirements:
1. **Input**: `xxx.pgm` (You can convert JPEG, GIF, and other formats to PGM ASCII format using IrfanView: [IrfanView](https://www.irfanview.com/))
2. **Output**: A compressed file `xxx.hc`
3. **Functions**:
    - **Compression**: Compress the PGM image using Huffman Coding.
      - You need to explain how to store the Huffman Tree and the associated codes in the file.
    - **Decompression**: Decompress the `xxx.hc` file back to its original format.
    - **Multiple File Compression**: Compress multiple files into a single file, supporting individual compression and decompression (optional bonus points).

---

#### Steps for Implementation:

1. **Compression**:
    - Read the input image file (`xxx.pgm` in ASCII format) and count the frequency of each pixel in the image.
    - Example command: 
      ```
      D:\> Huffman-Coding_Prog -c test.pgm
      ```
    - Create a Huffman coding tree using a Priority Queue based on the pixel frequencies.
    - Create a table of encodings for each pixel from the Huffman tree.
    - Encode the image using the generated codes and output the compressed image.
    
2. **Decompression**:
    - Read the encoded/compressed file (`xxx.hc`), decode it, and output the decompressed image.
    - Example command:
      ```
      D:\> Huffman-Coding_Prog -d test.hc testd.pgm
      ```

---

#### Required Output:
Your program must produce the following output:

1. **Pixel Frequency (Image Histogram)**: Display the frequency of each pixel value in the image.
2. **Pixel Value Encodings**: Provide a table listing each pixel value and its corresponding bit encoding, sorted from the smallest to the largest.
3. **File Sizes**: Display both the original image size and the compressed image size.
4. **Images**: Output both the compressed image (`xxx.hc`) and the uncompressed image (`testd.pgm`).
