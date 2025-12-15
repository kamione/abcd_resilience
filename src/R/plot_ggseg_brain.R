#' project values to a Desikan-Killiany brain atlas
#' @param dat A data frame.
#' @param fill A string.
#' @param title A string.
#' @param min A numeric.
#' @param max A numeric.
#' @returns A ggplot figure.

library(ggseg)
library(ggplot2)

plot_ggseg_brain <- function(
    dat, atlas = "dk", fill, min, max, break_int, title = NULL,
    fill_label = NULL, highlight_regions = NULL
) {
    colfunc <- colorRampPalette(
        c("#395D9C", "#358CA7", "white", "#F57A17", "#CE204E")
    )
    colors <- colfunc(100)
    
    if (!is.null(highlight_regions)) {
        dat <- dat |> 
            mutate(
                size = case_when(
                    label %in% highlight_regions ~ 1.5,
                    TRUE ~ 0.2
                ),
                color_group = case_when(
                    label %in% highlight_regions ~ "tomato3",
                    TRUE ~ "grey70"
                )
            )
        
        color_values <- dat %>% 
            pull(color_group) %>% 
            unique() %>% 
            sort()
    } else {
        dat <- dat %>% 
            mutate(
                size = 0.2, 
                color_group = "no"
            )
        color_values = "grey75"
    }

    if (atlas == "aseg") {
        fig <- dat |>
            ggplot() +
            geom_brain(
                atlas = aseg,
                position = position_brain("coronal"),
                mapping = aes(
                    fill = !!sym(fill),
                    color = color_group,
                    size = I(size)
                )
            ) +
            theme_void()
    } else if (atlas == "dk") {
        fig <- dk |>
            as_tibble() |>
            left_join(dat) |>
            ggplot() +
            geom_brain(
                atlas = dk,
                mapping = aes(
                    fill = !!sym(fill),
                    color = color_group,
                    size = I(size)
                )
            ) +
            theme_void()
    }
    fig <- fig +
        scale_fill_gradientn(
            colors = colors,
            na.value = "grey85",
            limits = c(min, max),
            breaks = seq(min, max, break_int)
        ) +
        scale_color_manual(values = color_values, guide = "none") +
        labs(title = title, fill = fill_label) +
        theme(
            legend.justification = c(0, 0.5)
        )

    return(fig)
}
