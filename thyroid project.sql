-- create database 
create database thyroiddb;

-- use created database
use thyroiddb

-- create a new table named 'probes'
create table probes (
id_ref varchar(20) primary key,
identifier varchar(100),
gene_symbol varchar(100),
gene_title varchar(225),
gene_id varchar(50),
chromosome_location varchar(50),
go_function varchar(50),
go_process varchar(1000),
go_component varchar(1000)
);

-- create second table named 'samnples'
create table samples (
sample_id varchar(20) primary key,
dose_kbg float,
absorbed_dose_gy float,
dose_group varchar(20)    -- control'/ 'low'/ 'medium'/ 'high'
);


-- create third table named 'expression'
create table expression ( 
id_ref varchar(20),
sample_id varchar(20),
expression_value float,
foreign key (id_ref) references probes(id_ref),
foreign key (sample_id) references samples(sample_id)
);

-- 1) Top 10 genes with the highest average expression overall
select p.gene_symbol, avg(e.expression_value) as avg_expr from expression e
join probes p on e.id_ref = p.id_ref
group by p.gene_symbol
order by avg_expr desc
limit 10;

-- 2) Genes with the biggest expression change between highest dose (17 Gy) and Control
select p.gene_symbol,
avg(case when s.dose_group = 'high' then e.expression_value
 end) as avg_high,
avg(case when s.dose_group ='control' then 
e.expression_value end) as avg_control,  
avg(case when s.dose_group = 'high' then e.expression_value 
end) - 
avg(case when s.dose_group = 'control'
 then e.expression_value end) as fold_diff 
from expression e 
join samples s on e.sample_id = s.sample_id
join probes p on e.id_ref = p.id_ref
group by p.gene_symbol 
order by ABS(fold_diff) desc
limit 20;

-- 3) Rank genes within each chromosome by average expression (window function)
select p.chromosome_location, p.gene_symbol,
avg(e.expression_value) as avg_expr,
rank() over (partition by p.chromosome_location order by
avg(e.expression_value) desc) as rnk
from expression e 
join probes p on e.id_ref = p.id_ref 
group by p.chromosome_location, p.gene_symbol;

-- 4) Expression trend of a gene across all 4 dose levels, ordered by dose
select p.gene_symbol, s.absorbed_dose_gy, e.expression_value
from expression e
join probes p on e.id_ref = p.id_ref
join samples s on e.sample_id = s.sample_id
where p.gene_symbol = 'agpat9'   -- swap in any gene of intrest
order by s.absorbed_dose_gy;

-- 5. TE to find genes with the largest expression range across dose levels (candidate dose-responsive genes)
with gene_range as (
select id_ref,
       max(expression_value) - min(expression_value) as
expr_range,
	avg(expression_value) as mean_expr
from expression
group by id_ref
)
select p.gene_symbol, gr.mean_expr, gr.expr_range
from gene_range gr
join probes p on gr.id_ref = p.id_ref
order by gr.expr_range desc
limit 25;

 -- 6. Genes linked to DNA damage / apoptosis response and their expression by dose group
 select p.gene_symbol, s.dose_group, avg(e.expression_value)as
 avg_expr
 from expression e
 join probes p on e.id_ref = p.id_REF	
 join samples s on e.sample_id = s.sample_id
 where p.go_process like '%dna damage%'or  p.go_process like 
 '%apoptotict%'
 group by p.gene_symbol, s.dose_group
 order by p.gene_symbol, s.dose_group;
 
 
 -- 7. Count of "Unknown" (unannotated) probes per chromosome
 select chromosome_location, count(*) as unknown_count
 from probes
 where gene_symbol ='unknown'
 group by chromosome_location
 order by unknown_count desc;