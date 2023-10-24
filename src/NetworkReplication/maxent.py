import pandas as pd
import numpy as np


def maxent(s_in, s_out, W):
    return np.outer(s_in, s_out.T) / W


def disaggregate(io: pd.core.frame.DataFrame, regional_sam: pd.core.frame.DataFrame):
    io_np = io.to_numpy()
    regional_sam_np = regional_sam.drop("Aggregate").to_numpy()

    arr = [list(io.index) + list(io.index), list(np.repeat(list(regional_sam.drop("Aggregate").index), len(io.index)))]
    idx = pd.MultiIndex.from_arrays(arr, names=["Sector", "Region"])
    df = pd.DataFrame(index=idx, columns=idx)
    for sector in io.index:
        for sector2 in io.index:
            W = io.loc[sector, sector2]
            s_out = W * regional_sam.drop("Aggregate").loc[:, sector]
            s_in = W * regional_sam.drop("Aggregate").loc[:, sector2]
            df.loc[sector, sector2] = maxent(s_in, s_out, W)
    return df


if __name__ == '__main__':
    d = {"Service": [500, 400, 100], "Agriculture": [200, 50, 150], "Manufacturing": [300, 150, 150]}
    regional_sam = pd.DataFrame(data=d, index=["Aggregate", "North", "South"])
    regional_sam.iloc[1:] = regional_sam.iloc[1:].div(regional_sam.iloc[0])
    io = pd.DataFrame(index=["Service", "Agriculture", "Manufacturing"],
                      data={"Service": [20, 30, 10], "Agriculture": [10, 80, 320], "Manufacturing": [60, 40, 200]})
    io = io / io.sum(axis=0)
    print(io, "\n")
    print(regional_sam)
    dio = disaggregate(io, regional_sam)
    df = pd.DataFrame(index=io.index, columns=io.index)
    for sector in io.index:
        for sector2 in io.index:
            print(sector, sector2)
            df.loc[sector, sector2] = dio.loc[sector, sector2].sum().sum()
