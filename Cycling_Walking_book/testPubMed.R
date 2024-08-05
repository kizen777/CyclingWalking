# 1. 必要なライブラリの読み込み
library(easyPubMed)
library(dplyr) 
library(stringr)
library(meta)
library(ggplot2)

# 2. PubMed APIキーの設定
set_entrez_key("a1b18c792b2447c1ded97f530279ecd72808")

# 3. 検索クエリの設定  
query <- "(cycling[Title/Abstract] OR bicycling[Title/Abstract] OR \"bicycle riding\"[Title/Abstract] OR walking[Title/Abstract]) AND (health[Title/Abstract] OR fitness[Title/Abstract] OR \"physical function\"[Title/Abstract] OR \"cardiovascular health\"[Title/Abstract] OR \"weight management\"[Title/Abstract] OR \"quality of life\"[Title/Abstract] OR \"mental health\"[Title/Abstract]) AND (comparison[Title/Abstract] OR comparative[Title/Abstract] OR versus[Title/Abstract] OR \"controlled trial\"[Publication Type] OR \"comparative study\"[Publication Type]) AND (\"2010\"[Date - Publication] : \"3000\"[Date - Publication])"

# 4. PubMed検索の実行
search_results <- entrez_search(db="pubmed", term=query, retmax=5000)
cat("検索された論文数:", length(search_results$ids), "\n")

# 5. データの取得
fetch_results <- entrez_fetch(db="pubmed", id=search_results$ids, rettype="abstract", parse=TRUE)

# 6. タイトルとアブストラクトの抽出
titles <- sapply(fetch_results, function(x) x$Title) 
abstracts <- sapply(fetch_results, function(x) x$AbstractText)
df <- data.frame(Title=titles, Abstract=abstracts, stringsAsFactors=FALSE)

# 7. 第1次スクリーニング
excl1 <- grepl("children|adolescents|review|meta-analysis|protocol", df$Title, ignore.case=TRUE)
after_excl1 <- df[!excl1, ]
cat("第1次スクリーニング後の論文数:", nrow(after_excl1), "\n") 

# 8. 第2次スクリーニング  
excl2 <- grepl("animal|in vitro|cell culture", after_excl1$Abstract, ignore.case=TRUE)
after_excl2 <- after_excl1[!excl2, ]  
cat("第2次スクリーニング後の論文数:", nrow(after_excl2), "\n")

# 9. Risk of Bias評価
# (ここではダミーのRisk of Bias評価を行っています。実際には各論文の詳細を確認して評価を行う必要があります）
set.seed(123)  
rob <- sample(c("Low", "Some concerns", "High"), nrow(after_excl2), replace=TRUE)
after_rob <- after_excl2[rob != "High", ]
cat("Risk of Bias評価後の論文数:", nrow(after_rob), "\n")

# 10. CONSORT図用データフレームの準備
dfMeta1 <- data.frame(
  Study = paste(after_rob$Title),
  Excl1 = ifelse(excl1[!excl2], "Excl1", ""), 
  Excl2 = ifelse(excl2, "Excl2", ""),
  RoB = ifelse(rob == "High", "RoB", "")
)

# 11. CONSORT図の作成
exclusions <- dfMeta1 %>% 
  summarise(
    Excl1 = sum(Excl1 != ""), 
    Excl2 = sum(Excl2 != ""),
    RoB = sum(RoB != "")
  )

flow_data <- data.frame(
  Stage = c("Studies identified through\ndatabase searching", 
            "Studies after duplicates removed",
            "Studies screened", 
            "Studies excluded",
            "Full-text articles\nassessed for eligibility", 
            "Full-text articles excluded,\nwith reasons",
            "Studies included in\nqualitative synthesis",
            "Studies included in\nquantitative synthesis\n(meta-analysis)"),
  Number = c(length(search_results$ids), 
             length(search_results$ids),
             length(search_results$ids),
             exclusions$Excl1,
             length(search_results$ids) - exclusions$Excl1,
             exclusions$Excl2 + exclusions$RoB,
             nrow(after_rob),
             nrow(after_rob))  
)

ggplot(flow_data, aes(x = 1, y = Number, label = Number)) +
  geom_text(nudge_x = 0.2, size = 3) +
  geom_segment(data = data.frame(y = c(length(search_results$ids) - exclusions$Excl1, 
                                       length(search_results$ids) - exclusions$Excl1 - exclusions$Excl2 - exclusions$RoB),
                                 yend = c(length(search_results$ids) - exclusions$Excl1 - exclusions$Excl2 - exclusions$RoB, nrow(after_rob)),
                                 x = 1, xend = 1),
               arrow = arrow(length = unit(0.2, "cm"))) +
  scale_y_continuous(breaks = flow_data$Number, labels = flow_data$Stage) +
  theme_void() +
  theme(axis.text.y = element_text(size = 10))

# 12. メタ分析の実施（ダミーデータを使用）
# (実際にはここで各研究のデータを抽出してメタ分析を行います）
set.seed(123)
effects <- rnorm(nrow(after_rob), 0.5, 0.2)  
meta_res <- metagen(effects, seq(0.1, 0.3, length.out=nrow(after_rob)))
forest(meta_res)

# 13. 結果の保存  
write.csv(dfMeta1, "CONSORT_data.csv", row.names=FALSE)