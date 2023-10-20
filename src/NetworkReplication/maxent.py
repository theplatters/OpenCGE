import pandas as pd
import numpy as np


def maxent(s_in, s_out, W):
    return np.outer(s_in, s_out.T) / W


def disaggregate(io: pd.core.frame.DataFrame, regional_sam: pd.core.frame.DataFrame):
    io_np = io.to_numpy()
    regional_sam_np = regional_sam.drop("Aggregate").to_numpy()
    s_in = np.reshape(np.sum(io_np, axis=0) * regional_sam_np, (-1, 1))

    print(s_in)
    s_out = np.reshape(np.sum(io_np, axis=1) * regional_sam_np, (-1, 1))
    print(s_out)
    return maxent(s_in, s_out,io.sum().sum())


if __name__ == '__main__':
    d = {"Service": [500, 400, 100], "Agriculture": [200, 50, 150], "Manufacturing": [300, 150, 150]}
    regional_sam = pd.DataFrame(data=d, index=["Aggregate", "North", "South"])
    regional_sam.iloc[1:] = regional_sam.iloc[1:].div(regional_sam.iloc[0])
    io = pd.DataFrame(index=["Service", "Agriculture", "Manufacturing"],
                      data={"Service": [20, 30, 10], "Agriculture": [10, 80, 320], "Manufacturing": [60, 40, 200]})
    io = io / io.sum(axis=0)
    print(io, "\n")
    print(regional_sam)
    print(disaggregate(io,regional_sam))
