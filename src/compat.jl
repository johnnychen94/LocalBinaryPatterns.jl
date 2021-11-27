if VERSION < v"1.5"
    # https://github.com/JuliaLang/julia/blob/70cc57cb36d839afe6ce56ea48ff6ed01bc262c4/base/int.jl#L524-L551
    function bitrotate(x::T, k::Integer) where {T <: Integer}
        (x << ((sizeof(T) << 3 - 1) & k)) | (x >>> ((sizeof(T) << 3 - 1) & -k))
    end
end
