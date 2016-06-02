#' Convert data.frame to cimg
#' @param mtrx a data.frame with (x, y, z) columns
#' @return cimg
#' @rdname morphology
#' @export
df2img = function(mtrx) {
    vars = c('x', 'y', 'z')
    mtrx = mtrx %>>% dplyr::select_(.dots=vars)
    .grid = dplyr::summarise_each_(mtrx, dplyr::funs(min, max), vars=vars) %>>%
        {expand.grid(x= seq(.$x_min, .$x_max),
                     y= seq(.$y_min, .$y_max),
                     z= seq(.$z_min, .$z_max))} %>>% dplyr::tbl_df()
    .grid %>>%
    dplyr::left_join(mtrx %>>% dplyr::mutate(v=1), by=vars) %>>%
    tidyr::replace_na(list(v=0)) %>>%
    reshape2::acast(x ~ y ~ z, dplyr::first, value.var='v', fill=0) %>>%
    {dim(.) = c(dim(.), 1); .} %>>%
    imager::as.cimg()
}

#' Convert cimg to data.frame
#' @param img a cimg object
#' @return a data.frame with (x, y, z) columns
#' @rdname morphology
#' @export
img2df = function(img) {
    as.array(img) %>>%
    {dim(.) = head(dim(.), 3); .} %>>%
    reshape2::melt(c('x', 'y', 'z')) %>>%
    as.data.frame %>>% dplyr::tbl_df()
}

#' Structuring element
#' @param coord a string
#' @param dimensions an integer
#' @return cimg structuring element
#' @rdname morphology
#' @export
get_se = function(coord=c('moore', 'neumann', 'hex'), dimensions=3) {
     coord = match.arg(coord)
     df = dplyr::data_frame(x=c(-1, 0, 1)) %>>%
         dplyr::mutate_(y=~x, z=~x) %>>%
         tidyr::expand_(dots=c('x', 'y', 'z'))
     if (coord == 'neumann') {
         df = dplyr::filter_(df, ~abs(x) + abs(y) + abs(z) < 2)
     } else if (coord == 'hex') {
         idx = trans_coord_hex(df) %>>% {.$x^2 + .$y^2 + .$z^2 < 1.1}
         df = dplyr::filter(df, idx)
     }
     if (dimensions < 3) {
         df = dplyr::filter_(df, ~z == 0)
     }
     df2img(df)
}

#' Filter img by surface
#' @param se a cimg object
#' @return cimg
#' @rdname morphology
#' @export
filter_surface = function(img, se) {
    img - imager::erode(img, se)
}

#' extract cells on suface
#' @return a logical vector
#' @rdname morphology
#' @export
detect_surface = function(mtrx, se) {
    axes = c('x', 'y', 'z')
    mins = dplyr::summarise_each_(mtrx, dplyr::funs(min), vars=axes)
    product = df2img(mtrx) %>>%
        filter_surface(se) %>>%
        img2df() %>>%
        dplyr::transmute_(
        x=~ x + mins$x - 1,
        y=~ y + mins$y - 1,
        z=~ z + mins$z - 1,
        surface=~ value > 0)
    dplyr::left_join(mtrx, product, by=axes)
}