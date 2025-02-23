
// 000 - 000, 0-0
// 001 - 100, 1-4
// 010 - 010, 2-2
// 011 - 110, 3-6
// 100 - 001, 4-1
// 101 - 101, 5-5
// 110 - 011, 6-3
// 111 - 111, 7-7


void WebRtcSpl_ComplexBitReverse(int16_t* __restrict complex_data, int stages) {
    int index = 0,                  // 当前的索引
        index_reversed = 0;         // 反转之后的索引
    int length = 1 << stages;
    int max = length - 1;

    // n 是 length
    // n 是 max number

    /* Decimation in time - re-order data */
    for (index = 1; index <= max; ++index) {
        int32_t* complex_data_ptr = (int32_t*)complex_data;
        int32_t temp = 0;

        /* Find out indexes that are bit-reversed. */
        int bit = length;                                                 // l = 8 = 1000, 它一定是 2^n，例如 0100, 0010, 0001
        do {
            bit >>= 1;                                                    // l = 4 = 0100
        } while (bit > max - index_reversed);                             // l > 7 - 4 = 3 = 0111 - 0100 = 0011
                                                                            // l > 7 - 2 = 5 = 0111 - 0010 = 0101
                                                                            // l > 7 - 1 = 6 = 0111 - 0001 = 0110
                                                                            // l > 7 - 0 = 7 = 0111 - 0000 = 0111
        // mask 的意义是位数，代表了第 N 位
        // index_reversed & (mask - 1) + mask，则是要把 > N 的位数改为 0，把第 N 位变成 1
        // 假设 mask = 0010，则 mask - 1 = 0001
        // 0101 & 0001 + 0010 = 0011
        int mask = bit - 1;
        index_reversed = (index_reversed & mask) + bit;

        if (index_reversed <= index) {
            continue;
        }

        /* Swap the elements with bit-reversed indexes.
            * This is similar to the loop in the stages == 7 or 8 cases.
            */
        temp = complex_data_ptr[index];  /* Real and imaginary */
        complex_data_ptr[index] = complex_data_ptr[index_reversed];
        complex_data_ptr[index_reversed] = temp;
    }
}