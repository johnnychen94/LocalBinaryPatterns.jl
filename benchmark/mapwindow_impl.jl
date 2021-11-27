using ImageFiltering

function lbp_origin_mapwindow(X)
    function _lbp(block)
        gc = block[2, 2]
        offsets = ((-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1))

        rst = 0
        for i in 1:length(offsets)
            o = offsets[i]
            rst += ifelse(gc <= block[2+o[1], 2+o[2]], 1, 0) << (i-1)
        end
        return rst
    end
    mapwindow(_lbp, X, (3, 3))
end
