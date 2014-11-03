
#define LH_DIVID_NUM 10
#define LH_DIVID_NUM_1 100
#define LH_DIVID_NUM_2 20
#define LH_MIN_THRESHOLD 16
#define LH_LEARN_RATE 0.001

#define LH_BRICK_STEP 4           // SubBrick's large size that used to compute
#define LH_PATCH_SIZE 1           // 提取Brick时的步宽，一定为1，别设大的，不然可能有问题，m_half_width_brick的计算与他有关。为1时可以保证没问题
#define LH_BRICK_PLY 3            // SubBrick's ply
#define LH_TARGET_DIM 5           // the target dimensions of onlinePCA subspace
#define LH_SAMPLE_SIZE 2          // Down-Sample rate

#define LH_T_MUST_UPDATE 200